# import numpy as np
# from JSB_tools.spe_reader import SPEFile
# from matplotlib import pyplot as plt
from openpyxl import load_workbook
import re

from pathlib import Path
shots_lis = []
shots_spe = []
for p in (Path.cwd()/'exp_data').iterdir():
    if p.is_dir():
        if p.name[-3:] == 'day':
            for f in p.iterdir():
                if m := re.match("shot([0-9]+)\.Lis", f.name):
                    shots_lis.append(int(m.groups()[0]))
                elif m := re.match("shot([0-9]+)\.Spe", f.name):
                    shots_spe.append(int(m.groups()[0]))


shots_lis = list(sorted(shots_lis))

for i in range(141):
    if i not in shots_lis:
        print(i)

print(shots_spe)
