import matplotlib.pyplot as plt
from JSB_tools.list_reader import MaestroListFile
import numpy as np
from lmfit.models import GaussianModel
from analysis import Shot
from scipy.signal import find_peaks
from uncertainties import unumpy as unp
from JSB_tools import Nuclide
from analysis import _get_maesto_list_shot_paths
from JSB_tools.list_reader import MaestroListFile
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
        print(shot)
        cal = EnergyCalMixin(erg_cal, load_erg_cal=False)
        cal.save_erg_cal(path=_get_maesto_list_shot_paths()[shot])

l = None
for shot in Shot.find_shots():
    if l is None:
        l = shot.list
    else:
        l += shot.list
        print()



