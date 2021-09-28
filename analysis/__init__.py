from __future__ import annotations

import warnings

import numpy as np
from JSB_tools.list_reader import MaestroListFile, get_merged_time_dependence
from pathlib import Path
import re
from openpyxl import load_workbook
from itertools import count
import pickle
from typing import Dict
from matplotlib import pyplot as plt
from functools import cache, cached_property
from JSB_tools import mpl_hist
from JSB_tools.spe_reader import SPEFile



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
        if isinstance(value, tuple) and len(value) == 2 and value[-1] is CONFIG_ATTRIB:
            self.__config_attribs__.append(key)
            self.__dict__[key] = value[0]
        else:
            self.__dict__[key] = value

    def __init__(self, shot_num):
        self.__config_attribs__ = []
        shot_metadata = ALL_SHOTS_METADATA[shot_num]

        self.flow = shot_metadata['flow'], CONFIG_ATTRIB
        self.he_flow = shot_metadata['He (SLPM)'], CONFIG_ATTRIB
        self.ar_flow = shot_metadata['Ar (SLPM)'], CONFIG_ATTRIB

        tube_len = shot_metadata['Length (Chamber-Filter m)']
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

        if self.shotnum in Shot.bad_shots:
            warnings.warn(f"Bad shot {self.shotnum} used!")

    @cached_property
    def list(self):
        try:
            path = _get_maesto_list_shot_paths()[self.shotnum]
        except KeyError:
            raise FileNotFoundError(f"No shot {self.shotnum}")
        out = MaestroListFile.from_pickle(path)
        time_offset(out)
        return out
 
    @staticmethod
    def find_shots(**attribs):
        shots = []
        for k, v in attribs.items():
            try:
                value = getattr(self, k)
            except:
                pass
            
        pass

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
            attribs: None to print all, 0 to print select information, or a list for custom information.

        Returns:

        """
        outs = []
        if attribs is None:
            attribs = self.__config_attribs__
        elif attribs == 0:
            attribs = ['flow', 'he_flow', 'ar_flow', 'foil_pos', 'cold_filter']

        for attrib in attribs:
            s = f'{attrib}={getattr(self, attrib)}'
            outs.append(s)
        return ', '.join(outs)



ALL_SHOTS_METADATA: Dict[int, dict] = dict()
ALL_SHOTS: Dict[int, Shot] = dict()


def __get_all_shots(load=False):
    global ALL_SHOTS_METADATA
    if load:
        with open('excel.pickle', 'rb') as f:
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
    with open('excel.pickle', 'wb') as f:
        pickle.dump(ALL_SHOTS_METADATA, f)


__get_all_shots(True)

if __name__ == '__main__':
    s = SPEFile('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/Nickel/Nickel.Spe')
    s.plot_erg_spectrum()
    shots = list(map(lambda x: Shot(x).list, [26, 27]))
    ax = None
    sig, bg, bins = get_merged_time_dependence([Shot(26).list, Shot(27).list], 218.5)
    mpl_hist(bins, sig)
    ax = Shot(134).list.plot_erg_spectrum(time_max=60*4, remove_baseline=True)
    Shot(34).list.plot_erg_spectrum(time_max=60*4, ax=ax, remove_baseline=True)

    plt.show()



