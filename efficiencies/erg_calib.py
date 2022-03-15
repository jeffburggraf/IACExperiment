import matplotlib.pyplot as plt
from JSB_tools.maestro_reader import MaestroListFile
import numpy as np
from lmfit.models import GaussianModel
from analysis import Shot
from scipy.signal import find_peaks
from uncertainties import unumpy as unp
from JSB_tools import Nuclide
from analysis import _get_maesto_list_shot_paths
from JSB_tools.maestro_reader import MaestroListFile
from JSB_tools.spe_reader import SPEFile
from pathlib import Path
import re
from JSB_tools.spe_reader import EnergyCalMixin


cwd = Path(__file__).parent

class_dir = cwd/'NewErgCals'


def get_erg_cal(p):
    with open(p) as f:
        lines = f.readlines()
    return list(map(float, lines[lines.index('$ENER_FIT:\n') + 1].split()))

d = {}
if True:
    for path in class_dir.iterdir():
        if m := re.match("Shot ([0-9]+)-([0-9]+)-", path.name):
            shots = list(range(int(m.groups()[0]), int(m.groups()[1]) + 1))
        elif m := re.match("Shot ([0-9]+)-", path.name):
            shots = [int(m.groups()[0])]
        elif re.match("group", path.name):
            with open(path) as f:
                line = f.readlines()[1]
            m = re.match(' ?Shots (.+)', line)
            assert m
            shots = eval(m.groups()[0])
        else:
            continue

        erg_cal = get_erg_cal(path)

        for shot in shots:
            d[shot] = erg_cal
print(d)

        # for shot in shots:
        #     print(shot)
        #     cal = EnergyCalMixin(erg_cal, load_erg_cal=False)
        #     cal.save_erg_cal(path=_get_maesto_list_shot_paths()[shot])

#
# #  =========================================
# imax = 100
# min_time = 350
# #  =========================================
# i = 0
# l = None
# for shot in Shot.find_shots():
#
#     if shot.max_time < min_time:
#         continue
#     if l is None:
#         l = shot.list
#     else:
#         l.__iadd__(shot.list, fast_calc=True)
#     print(shot)
#     i += 1
#     if i > imax:
#         break
#
# l.plotly(erg_min=50, erg_max=2400, time_bin_width=7.5, time_step=3,convolve_overlay_sigma=1)
#


