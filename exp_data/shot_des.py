import re
from pathlib import Path
from openpyxl import load_workbook
#
wb = load_workbook(filename = Path(__file__).parent/'IAC Run Spreadsheet.xlsx')
sheet = wb['Sheet1']


class ShotInfo:
    def __init__(self, entries_dict):
        self.entries_dict = entries_dict
        self.shot_num = self.entries_dict['Run #']
        self.n_pulses = self.entries_dict['# pulses']
        if self.n_pulses is None:
            self.is_valid = False
        else:
            self.is_valid = True

        flow = [''.join(map(str,
                            (entries_dict[f'Upstream{i}'], entries_dict[f'Center{i}'],
                             entries_dict[f'Downstream{i}'])))
                for i in [1, 2]]
        self.flow = flow
        self.entries_dict['flow'] = self.flow


all_shots = {}
header = []
n = None
for h in sheet['2']:
    h = h.value
    if not h:
        break
    if 'Upstream' in h:
        if n is None:
            n = 1
        else:
            n += 1
    if re.match("Upstream|Downstream|Center", h):
        header.append(f'{h}{n}')
    else:
        header.append(h)

print(header)

for i in range(3, 210):
    entries = [c.value for c in sheet[f'{i}']]
    info = {}
    for e, h in zip(entries, header):

        info[h] = e
    try:
        all_shots[info['Run #']] = ShotInfo(info)
    except KeyError:
        pass

