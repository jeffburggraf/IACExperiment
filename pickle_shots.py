from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np

from JSB_tools.list_reader import MaestroListFile


data_dir = Path(__file__).parent/'exp_data'
# data_dir = Path('/Volumes/NO NAME')


def get_maesto_list_shot_paths():
    out = {}
    for path in data_dir.iterdir():
        if path.is_dir() and re.match('.+day', path.name):
            for path in path.iterdir():
                if m := re.match(r'shot([0-9]+)\.Lis', path.name):
                    out[int(m.groups()[0])] = path
    out = {k:v for k, v in sorted(out.items(), key=lambda x: x[0])}

    return out


def get_mpant_mca_shot_paths():
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





if __name__ == '__main__':
    # p132 = get_maesto_list_shot_paths()[132]
    # l = MaestroListFile(p132)
    # l.SPE.plot_erg_spectrum()
    # plt.show()
    times, files = [], []
    for path in get_maesto_list_shot_paths().values():
        # m = MaestroListFile(path)
        m = MaestroListFile.from_pickle(path)
        t = m.total_real_time
        i = np.searchsorted(-np.array(times), -t)
        times.insert(i, t)
        files.insert(i, path.name)
        print(f'Longest: {files[:15]}\n {times[:15]} s\n')
        # m.pickle()
