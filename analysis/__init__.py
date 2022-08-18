from __future__ import annotations
import warnings
import numpy as np
from JSB_tools.maestro_reader import MaestroListFile, get_merged_time_dependence
from pathlib import Path
import re
from openpyxl import load_workbook
from itertools import count
import pickle
from typing import Dict, List, Union
from matplotlib import pyplot as plt
from functools import cache, cached_property
from JSB_tools import mpl_hist
from JSB_tools.spe_reader import SPEFile
from mpant_reader import MPA
from datetime import datetime


cwd = Path(__file__).parent


def time_offset(l: MaestroListFile):
    times = l.realtimes[np.where(l.sample_ready_state == 1)]
    time = np.median(times)
    l.time_offset(-time)


data_dir = Path(__file__).parent.parent/'exp_data'


class Default:
    pass

_FF_EXP_DECAY_RATES = None

def get_yield(times):
    pass


@cache
def _get_maesto_list_shot_paths():
    out = {}
    for path in data_dir.iterdir():
        if path.is_dir() and re.match('.+day', path.name):
            for path in path.iterdir():
                if m := re.match(r'shot([0-9]+)\.Lis', path.name):
                    out[int(m.groups()[0])] = path
    out = {k:v for k, v in sorted(out.items(), key=lambda x: x[0])}

    return out

@cache
def _get_mpant_mca_shot_paths():
    out = {}
    for path in data_dir.iterdir():
        if path.is_dir() and re.match('.+day', path.name):
            mca_path = path/'MCA'
            if not mca_path.exists():
                continue
            for path in mca_path.iterdir():
                if m := re.match(r'shot([0-9]+)\.mpa', path.name):
                    out[int(m.groups()[0])] = path
    out = {k:v for k, v in sorted(out.items(), key=lambda x: x[0])}
    return out


# CONFIG_ATTRIB is used as an indicator to set which Shot attributes will be used in __eq__ to test if two run
# configurations have the same configuration.
CONFIG_ATTRIB = object


class Shot:
    bad_shots = [6, 42, 43, 44, 82, 119, 123, 136, 140]
    # equal_test_attribs = {'flow', 'mylar', }

    shape_cal = [0.4448941791251471, 0.000246982583338865, -1.538593499008935e-08]

    @staticmethod
    def nickel_list():
        p = data_dir/'Nickel'/'Nickel.Lis'
        try:
            return MaestroListFile.from_pickle(p)
        except FileNotFoundError:
            return MaestroListFile(p)

    @staticmethod
    def nickel_spe():
        p = data_dir / 'Nickel' / 'Nickel.Spe'
        return SPEFile(p)

    def __setattr__(self, key, value):
        """
        If attribute is set with a Tuple for which the second value is CONFIG_ATTRIB, then this attribute is considered
            in determining if two shots are equal. These attributes are also the default used in __repr__
        Args:
            key:
            value:

        Returns:

        """
        if isinstance(value, tuple) and len(value) == 2 and value[-1] is CONFIG_ATTRIB:
            self.__config_attribs__.append(key)
            self.__dict__[key] = value[0]
        else:
            self.__dict__[key] = value

    @staticmethod
    def background_spe(llnl=True):
        if llnl:
            return SPEFile(cwd.parent/'exp_data'/'tuesday'/'BG.Spe')
        return MPA(cwd.parent/'exp_data'/'tuesday'/"MCA"/'BG001.Spe')

    def __init__(self, shot_num, load_erg_cal=True):
        """

        """
        self.__config_attribs__ = []
        shot_metadata = ALL_SHOTS_METADATA[shot_num]
        self.load_erg_cal = load_erg_cal

        self.flow = shot_metadata['flow'], CONFIG_ATTRIB
        self.he_flow = shot_metadata['He (SLPM)'], CONFIG_ATTRIB
        self.ar_flow = shot_metadata['Ar (SLPM)'], CONFIG_ATTRIB

        self.converter = shot_metadata['Converter'], CONFIG_ATTRIB

        try:
            self.pressure = shot_metadata['P (psai)']/14.5
        except TypeError:
            self.pressure = 0

        tube_len = shot_metadata['Tube length (m)']
        try:
            self.tube_len = eval(tube_len.replace('=', '')), CONFIG_ATTRIB
        except AttributeError:
            self.tube_len = float(tube_len), CONFIG_ATTRIB

        self.cold_filter = shot_metadata['Filter temp'] == 'LN2', CONFIG_ATTRIB
        self.mylar = shot_metadata['Mylar (um)'], CONFIG_ATTRIB

        foil_pos = shot_metadata['#1 pos from upstream']
        self.foil_pos = "center" if foil_pos == 'center' else 'upstream', CONFIG_ATTRIB

        self.flow_stop = shot_metadata['Flow stop (s)'], CONFIG_ATTRIB

        if self.flow_stop is None:
            self.flow_stop = 0

        self.num_filters = shot_metadata['#'], CONFIG_ATTRIB
        self.beam_duration = shot_metadata['Duration (s)'], CONFIG_ATTRIB

        self.n_pulses = int(shot_metadata['# pulses'])
        self.comment = shot_metadata['Comments']
        self.shotnum = shot_metadata['Run #']

        # for k, v in shot_metadata.items():
        #     print(k.__repr__(), v)
        self.iac_filter_pos = shot_metadata['IAC Det Pos.']
        self.llnl_filter_pos = shot_metadata['LLNL Det Pos.'], CONFIG_ATTRIB

        if self.shotnum in Shot.bad_shots:
            warnings.warn(f"Bad shot {self.shotnum} used!")
        self.good_shot = self.shotnum not in Shot.bad_shots

    @staticmethod
    def get_sigma(erg):
        return np.sum([c*erg**i for i, c in enumerate(Shot.shape_cal)])

    @staticmethod
    def get_fwhm(erg):
        return 2.355*Shot.get_sigma(erg)

    @staticmethod
    def get_peak_width(erg):
        return 4*Shot.get_sigma(erg)

    @property
    def max_time(self):
        return self.list.times[-1]

    @property
    def filters_before_llnl_det(self):
        return self.llnl_filter_pos - 1

    @property
    def filters_before_iac_det(self):
        return self.iac_filter_pos - 1

    @cached_property
    def list(self) -> MaestroListFile:
        try:
            path = _get_maesto_list_shot_paths()[self.shotnum]
        except KeyError:
            raise FileNotFoundError(f"No shot {self.shotnum}")
        try:
            out = MaestroListFile.from_pickle(path, load_erg_cal=self.load_erg_cal)
        except FileNotFoundError:
            out = MaestroListFile(path)
            out.pickle()
        time_offset(out)
        out.allow_pickle = False
        return out

    @property
    def llnl_spe(self):
        return self.list.SPE

    @cached_property
    def iac_spe(self) -> SPEFile:
        try:
            path = _get_mpant_mca_shot_paths()[self.shotnum]
        except KeyError:
            raise FileNotFoundError(f"No shot {self.shotnum}.mpa")
        try:
            out = MPA.from_pickle(path)
        except FileNotFoundError:
            out = MPA(path)
        return out

    valid_values = {'foil_pos': ['center', 'upstream'],
                    'mylar': [0, 2.5, 5, 10, 20],
                    'flow_stop': [0, 10, 20, 40],
                    'tube_len': [4.16, 6.14, 9.44, 12.64],
                    'flow': ['100001', '111111', '010010', '111010'],
                    'beam_duration': [3, 5, 20],
                    'num_filters': [2, 3, 4],
                    'llnl_filter_pos': [2, 1]}

    @staticmethod
    def find_shots(good_shot_only=True, eval_func='1==1', flow=None, he_flow=None, ar_flow=None, tube_len=None,
                   cold_filter=None, mylar=None, foil_pos=None, flow_stop=None, num_filters=None, beam_duration=None,
                   converter='W', llnl_filter_pos=1) -> List[Shot]:
        """
        None always means not to use as search criteria.

        Args:
            good_shot_only: True means dont include "bad" shots.

            eval_func: string of python code to be evaluated with 'self' bound to current instance.
                e.g.: self.he_flow + self.ar_flow >= 1.0

            flow: Flow patter. Options: 010010 (default), 100001 (offset), 111111 (All in all out),
                or 111010 (All in one out)

            he_flow: He flow in SLPM
            ar_flow: Ar flow in SLPM

            tube_len: None for any. Options are: 4.16, 6.14, 9.44, 12.64  (in meters)
                Can be Tuple of len 1, e.g. (1,) to refer to length of the second tube length, as in, [4.16, 6.14, 9.44, 12.64][1]

            cold_filter: True or False
            mylar: 0 for no mylar. options are: 0, 2.5, 5, 7.5, 10, 20
            foil_pos: "1cm" for upstream. , "center" for center.
            flow_stop: 0 for no stopping of gas. Otherwise, number of seconds before stop.
            num_filters: Number of filters used.
            beam_duration: Seconds beam was running. Usually 3.

            converter: W for tungsten converter.

        Returns:

        """
        # shots = []
        if isinstance(tube_len, tuple) and len(tube_len) == 1 and type(tube_len[0]) == int:
            tube_len = [4.16, 6.14, 9.44, 12.64][tube_len[0]]

        attribs = {'flow': flow, 'he_flow': he_flow, 'ar_flow': ar_flow, 'tube_len': tube_len,
                   'cold_filter': cold_filter, 'mylar': mylar, 'foil_pos': foil_pos, 'flow_stop': flow_stop,
                   'num_filters': num_filters, 'beam_duration': beam_duration, 'converter': converter,
                   'llnl_filter_pos': llnl_filter_pos}

        for k, v in attribs.items():
            if v is None:
                continue
            try:
                assert v in Shot.valid_values[k], f'Invalid value, "{v}", for argument {k}'
            except KeyError:
                continue

        attribs = dict(filter(lambda k_v: k_v[1] is not None, attribs.items()))

        if eval_func != '1==1':
            assert isinstance(eval_func, str), "`eval_func` must be a string"
            assert 'self' in eval_func, "Bad eval_func. Example usage: 'self.ar_flow + self.he_flow >= 1.0'"

        eval_func = eval_func.replace('self', 'shot')

        for shot_num in ALL_SHOTS_METADATA:

            if good_shot_only and shot_num in Shot.bad_shots:
                continue
            shot = Shot(shot_num)
            try:
                if all([getattr(shot, name) == value for name, value in attribs.items()]):

                    if eval(eval_func):
                        # shots.append(shot)
                        yield shot

            except AttributeError:
                raise

        # return shots

    @staticmethod
    def sum_shots_list(shots_list: List[Union[Shot, int]], truncate=False) -> MaestroListFile:
        out = None

        for shot in shots_list:
            if isinstance(shot, int):
                shot = Shot(shot)
            l = shot.list
            if out is None:
                out = l
            else:
                out.__iadd__(l, truncate)
        return out

    @property
    def acq_time(self):
        with open(Path(__file__).parent/'shot_times.pickle', 'rb') as f:
            times = pickle.load(f)

        return times[self.shotnum]

    @staticmethod
    def n_shots_array(shots_list: List[Union[Shot, int]], time_bins):
        """
        Return the number of shots in list `shot_list` which are acquiring data for each time interval in time_bins
        Args:
            shots_list:
            time_bins:

        Returns:

        """
        def f(_tmax, _time_bins):
            _out = np.zeros(len(_time_bins) - 1, dtype=float)
            if _tmax >= _time_bins[-1]:
                return np.ones_like(_out)
            if _time_bins[0] >= _tmax:
                return _out
            n = np.searchsorted(_time_bins, _tmax, side='right') - 1
            _out[:n] += 1
            _out[n] += (_tmax - _time_bins[n])/(_time_bins[n+1] - _time_bins[n])
            return _out

        out = None
        for shot in shots_list:
            if isinstance(shot, int):
                shot = Shot(shot)
            r = f(shot.acq_time, time_bins)
            if out is None:
                out = r
            else:
                out += r

        return out

    def __eq__(self, other: Shot):
        assert isinstance(other, Shot)
        for key in self.__config_attribs__:
            v1 = getattr(self, key)
            v2 = getattr(other, key)
            if v1 != v2:
                return False
        return True

    def __repr__(self, attribs='all'):
        """

        Args:
            attribs: 'all' to print all shot meta data, 0 to print select information, or a list for custom attributes.

        Returns:

        """
        outs = []
        if attribs == 'all':
            attribs = ['shotnum'] + self.__config_attribs__ +  ['comment']
        elif attribs == 0:
            attribs = ['shotnum', 'flow', 'he_flow', 'ar_flow', 'foil_pos', 'cold_filter', 'comment']
        else:
            if hasattr(attribs, '__iter__') and all(isinstance(x, str) for x in attribs):
                pass
            else:
                raise ValueError(f'Invalid `attribs` argument,{attribs}')

        for attrib in attribs:
            s = f'{attrib}={getattr(self, attrib)}'
            outs.append(s)
        return ', '.join(outs)


ALL_SHOTS_METADATA: Dict[int, dict] = dict()
# ALL_SHOTS: Dict[int, Shot] = dict()


def __get_all_shots_data(load=False):
    global ALL_SHOTS_METADATA
    if load:
        with open(cwd/'excel.pickle', 'rb') as f:
            ALL_SHOTS_METADATA = pickle.load(f)
            return ALL_SHOTS_METADATA

    # Make header list
    wb = load_workbook(filename=data_dir / 'IAC Run Spreadsheet.xlsx', data_only=True)
    sheet = wb['Sheet1']
    header_list = []
    for h in sheet['2']:
        try:
            header_list.append(h.value.rstrip().lstrip())
        except AttributeError:
            continue
    header_list.append("")

    comment = None  # Keep comment persistent since it only occurs once in merged rows.
    for i in count(3):  # data starts in row 3.
        shot_metadata = {h: v.value for h, v in zip(header_list, sheet[f'{i}'])}
        if shot_metadata['Comments'] is not None:
            comment = shot_metadata['Comments']
        shot_metadata['Comments'] = comment

        flow = ''
        for k in ['Upstream Inlet', 'Center Inlet', 'Downstream Inlet', 'Upstream Outlet', 'Center Outlet',
                  'Downstream Outlet']:
            flow += str(shot_metadata.pop(k))

        shot_metadata['flow'] = flow

        run_num = shot_metadata.get('Run #', None)
        if run_num is None:  # end of data
            break
        run_num = int(run_num)
        ALL_SHOTS_METADATA[run_num] = shot_metadata
    ALL_SHOTS_METADATA = {shot_num: {k.rstrip().lstrip(): v for k, v in meta.items()}  for shot_num, meta in ALL_SHOTS_METADATA.items()}

    with open(cwd/'excel.pickle', 'wb') as f:
        pickle.dump(ALL_SHOTS_METADATA, f)


__get_all_shots_data(False)


def save_shot_times():
    shots_times = {}
    for shot in Shot.find_shots():
        shots_times[shot.shotnum] = shot.list.times[-1]
        #
    shots_times = {k: v for k, v in sorted(shots_times.items(), key=lambda k_v: -k_v[-1])}
    with open(Path(__file__).parent/'shot_times.pickle', 'wb') as f:
        pickle.dump(shots_times, f)


if __name__ == '__main__':
    # save_shot_times()
    shot = Shot.n_shots_array(None, None, )
    # s = Shot(134, )
    # s.list.plotly()



