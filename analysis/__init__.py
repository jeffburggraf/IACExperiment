from __future__ import annotations

import warnings

import numpy as np
from JSB_tools.list_reader import MaestroListFile, get_merged_time_dependence
from pathlib import Path
import re
from openpyxl import load_workbook
from itertools import count
import pickle
from typing import Dict, List
from matplotlib import pyplot as plt
from functools import cache, cached_property
from JSB_tools import mpl_hist
from JSB_tools.spe_reader import SPEFile
from mpant_reader import MPA

cwd = Path(__file__).parent


def time_offset(l: MaestroListFile):
    times = l.realtimes[np.where(l.sample_ready_state == 1)]
    time = np.median(times)
    l.time_offset(-time)


data_dir = Path(__file__).parent.parent/'exp_data'


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
    bad_shots = [42, 43, 44, 82, 123, 136]
    equal_test_attribs = {'flow', 'mylar', }

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
    def background_spe():
        return SPEFile(cwd.parent/'exp_data'/'tuesday'/'BG.Spe')

    def __init__(self, shot_num):
        """

        """
        self.__config_attribs__ = []
        shot_metadata = ALL_SHOTS_METADATA[shot_num]

        self.flow = shot_metadata['flow'], CONFIG_ATTRIB
        self.he_flow = shot_metadata['He (SLPM)'], CONFIG_ATTRIB
        self.ar_flow = shot_metadata['Ar (SLPM)'], CONFIG_ATTRIB

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
        self.num_filters = shot_metadata['#'], CONFIG_ATTRIB
        self.beam_duration = shot_metadata['Duration (s)'], CONFIG_ATTRIB
        self.convertor = shot_metadata['Converter'], CONFIG_ATTRIB

        self.n_pulses = int(shot_metadata['# pulses'])
        self.comment = shot_metadata['Comments']
        self.shotnum = shot_metadata['Run #']

        # for k, v in shot_metadata.items():
        #     print(k.__repr__(), v)
        self.iac_filter_pos = shot_metadata['IAC Det Pos.']
        self.llnl_filter_pos = shot_metadata['LLNL Det Pos.']

        if self.shotnum in Shot.bad_shots:
            warnings.warn(f"Bad shot {self.shotnum} used!")

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
    def list(self):
        try:
            path = _get_maesto_list_shot_paths()[self.shotnum]
        except KeyError:
            raise FileNotFoundError(f"No shot {self.shotnum}")
        try:
            out = MaestroListFile.from_pickle(path)
        except FileNotFoundError:
            out = MaestroListFile(path)
            out.pickle()
        time_offset(out)
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
        return MPA(path)

    @staticmethod
    def find_shots(**attribs) -> List[Shot]:
        shots = []
        for shot_num in ALL_SHOTS_METADATA:
            if shot_num in Shot.bad_shots:
                continue
            shot = Shot(shot_num)
            try:
                if all([getattr(shot, name) == value for name, value in attribs.items()]):
                    shots.append(shot)
            except AttributeError:
                raise
        return shots

    def __eq__(self, other: Shot):
        assert isinstance(other, Shot)
        for key in self.__config_attribs__:
            v1 = getattr(self, key)
            v2 = getattr(other, key)
            if v1 != v2:
                return False
        return True

    def __repr__(self, attribs=0):
        """

        Args:
            attribs: 'all' to print all shot meta data, 0 to print select information, or a list for custom attributes.

        Returns:

        """
        outs = []
        if attribs == 'all':
            attribs = self.__config_attribs__
        elif attribs == 0:
            attribs = ['shotnum', 'comment', 'flow', 'he_flow', 'ar_flow', 'foil_pos', 'cold_filter']
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
    wb = load_workbook(filename=data_dir / 'IAC Run Spreadsheet.xlsx')
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

        flow = ''
        for k in ['Upstream Inlet', 'Center Inlet', 'Downstream Inlet', 'Upstream Outlet', 'Center Outlet', 'Downstream Outlet']:
            flow += str(shot_metadata.pop(k))

        shot_metadata['flow'] = flow
        shot_metadata['Comments'] = comment

        run_num = shot_metadata.get('Run #', None)
        if run_num is None:  # end of data
            break
        run_num = int(run_num)
        ALL_SHOTS_METADATA[run_num] = shot_metadata
    ALL_SHOTS_METADATA = {shot_num: {k.rstrip().lstrip(): v for k, v in meta.items()}  for shot_num, meta in ALL_SHOTS_METADATA.items()}

    with open(cwd/'excel.pickle', 'wb') as f:
        pickle.dump(ALL_SHOTS_METADATA, f)


__get_all_shots_data(False)

if __name__ == '__main__':
    for s in Shot.find_shots(flow='010010', tube_len=12.64, num_filters=2, cold_filter=False, he_flow=0.5, ar_flow=0.5,
                             mylar=0, foil_pos='upstream'):
        print(s)
    # ax = Shot(56).list.plot_time_dependence(218.5)
    # Shot(57).list.plot_time_dependence(218.5, ax=ax)
    # Shot(58).list.plot_time_dependence(218.5, ax=ax)
    # plt.show()
    # Shot().tube_len
    # s = SPEFile('/exp_data/Nickel/Nickel_original.Spe')
    # s.plot_erg_spectrum()
    # shots = list(map(lambda x: Shot(x).list, [26, 27]))
    # ax = None
    # sig, bg, bins = get_merged_time_dependence([Shot(26).list, Shot(27).list], 218.5)
    # mpl_hist(bins, sig)
    # ax = Shot(134).list.plot_erg_spectrum(time_max=60*4, remove_baseline=True)
    # Shot(34).list.plot_erg_spectrum(time_max=60*4, ax=ax, remove_baseline=True)
    #
    # plt.show()



