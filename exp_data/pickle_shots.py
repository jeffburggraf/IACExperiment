"""
Should do it all except for efficiency
"""
from analysis import _get_maesto_list_shot_paths
from pathlib import Path
import re
from multiprocessing import Process
from JSB_tools.maestro_reader import MaestroListFile
from JSB_tools.spe_reader import EnergyCalMixin

cwd = Path(__file__).parent

class_dir = cwd / 'NewErgCals'


def get_erg_cal(p):
    with open(p) as f:
        lines = f.readlines()
    return list(map(float, lines[lines.index('$ENER_FIT:\n') + 1].split()))


erg_cals = {}

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
            erg_cals[shot] = erg_cal


def pickle_list(shots_paths):
    for shot, path in shots_paths:
        if not path.exists():
            continue# _get_maesto_list_shot_paths().items():
        s = MaestroListFile(path)
        try:
            s.erg_calibration = erg_cals[shot]
        except KeyError:
            pass
        s.pickle(path)
        print(f"Shot {shot}")


if __name__ == '__main__':
    l = list(_get_maesto_list_shot_paths().items())
    #  ========================================================
    multiprocessess = True
    min_erg=40
    max_erg= 2700
    #  ========================================================

    if multiprocessess:
        n_pools = 6

        di = len(l)//6
        max_i = max(_get_maesto_list_shot_paths().keys())
        processes = []
        pnickel= Process()
        for i1, i2 in zip(range(0, len(l), di), range(di, len(l)+di, di)):
            # i2 = max()
            p1 = Process(target=pickle_list, args=(l[i1: i2],))
            processes.append(p1)

        for p in processes:
            p.start()

        for p in processes:
            p.join()

        print("done!")

    else:
        pickle_list(l)

