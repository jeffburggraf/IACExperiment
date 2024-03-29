import re
from pathlib import Path
from openpyxl import load_workbook
from typing import Dict, List
from typeguard import check_type

#
wb = load_workbook(filename=Path(__file__).parent/'exp_data'/'IAC Run Spreadsheet.xlsx')
sheet = wb['Sheet1']


def check_attribs(attrib, value):
    try:
        t = {'flow_pat': List[str]}[attrib]
    except KeyError:
        return
    check_type(attrib, value, t)  #  raises error


def get_all_shots_with(**attribs):
    def test(name, value, shot_metadata: ShotInfo):
        check_attribs(name, value)
        try:
            other = getattr(shot_metadata, name)
        except AttributeError:
            try:
                other = shot_metadata.entries_dict[name]
            except KeyError as e:
                print(e)
                return False
        if other == value:
            return True

        return False
    outs = []
    for shot in all_shots.values():
        if all(test(name, value, shot) for name, value in attribs.items()):
            outs.append(shot)
    return outs


class ShotInfo:
    def __init__(self, entries_dict):
        self.entries_dict = entries_dict
        self.shot_num = self.entries_dict['Run #']
        self.n_pulses = self.entries_dict['# pulses']
        self.he_flow = entries_dict['He (SLPM)']
        self.ar_flow = entries_dict['Ar (SLPM)']
        self.mylar_thickness = entries_dict['Mylar (um)']
        self.ufoil_pos = entries_dict['foil pos']

        if self.n_pulses is None:
            self.is_valid = False
        else:
            self.is_valid = True

        flow = ''.join(map(str,
                        (entries_dict[f'Upstream Inlet'], entries_dict['Center Inlet'], entries_dict['Downstream Inlet']
                         ,entries_dict[f'Upstream Outlet'], entries_dict['Center Outlet'],
                         entries_dict['Downstream Outlet'])))

        self.flow = flow
        self.entries_dict['flow_pat'] = self.flow
        # self.entries_dict['foil pos'] = self.ufoil_pos
        if self.entries_dict['Filter temp'] is None:
            self.entries_dict['Filter temp'] = 'Room'
        self.is_ln2 = True if self.entries_dict['Filter temp'] == 'LN2' else False

    def __repr__(self, labels=None):
        if labels is None:
            labels = ['He (SLPM)', 'Ar (SLPM)', 'flow_pat', 'Mylar (um)', 'Filter temp', 'foil pos']
        outs = [f"Shot: {self.shot_num}"]
        for k, v in self.entries_dict.items():
            if k in labels:
                outs.append(f"{k}: {v}")

        return '; '.join(outs)


header = []

for h in sheet['2']:
    h = h.value
    if not h:
        break
    header.append(h)


print(header)
all_shots: Dict[int, ShotInfo] = {}

for i in range(3, 210):
    entries = [c.value for c in sheet[f'{i}']]
    info = {}
    for e, h in zip(entries, header):
        info[h] = e
    try:
        info['foil pos'] = info['#1 pos from upstream']
        info.pop('#1 pos from upstream')
        all_shots[info['Run #']] = ShotInfo(info)
    except KeyError:
        raise

if __name__ == '__main__':
    for i in get_all_shots_with(flow=['100', '001']):
        print(i)