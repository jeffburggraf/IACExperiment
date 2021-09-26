import numpy as np
from JSB_tools.list_reader import MaestroListFile
from pathlib import Path
import re
from openpyxl import load_workbook
from itertools import count
import pickle
from typing import Dict
from matplotlib import pyplot as plt
from functools import cache, cached_property
from JSB_tools import mpl_hist


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


class Shot:
    bad_shots = [42, 43, 44, 82, 136]

    def __init__(self, shot_num):
        shot_metadata = ALL_SHOTS_METADATA[shot_num]
        self.shotnum = shot_metadata['Run #']
        self.flow = shot_metadata['flow']
        self.mylar = shot_metadata['Mylar (um)']
        self.n_pulses = int(shot_metadata['# pulses'])
        self.comment = shot_metadata['Comments']

    @cached_property
    def list(self):
        path = _get_maesto_list_shot_paths()[self.shotnum]
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


ALL_SHOTS_METADATA: Dict[int, dict] = dict()
ALL_SHOTS: Dict[int, Shot] = dict()


def get_all_shots(load=False):
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

    with open('excel.pickle', 'wb') as f:
        pickle.dump(ALL_SHOTS_METADATA, f)


get_all_shots(True)

if __name__ == '__main__':
    import timeit
    shot = Shot(42)
    # ax = shot.list.plotly()
    shot.list.plot_erg_spectrum()
    s , _, bins = shot.list.get_time_dependence(218)
    ax = mpl_hist(bins, s)
    # Shot(139).list.plot_percent_live()
    # Shot(136).list.plot_percent_live()
    # Shot(136).list.plot_count_rate()
    # Shot(122).list.plot_percent_live()
    # Shot(122).list.plot_erg_spectrum(ax=ax)
    plt.show()
    # s = Shot(132)



