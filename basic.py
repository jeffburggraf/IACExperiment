import warnings

import numpy as np

from JSB_tools.list_reader import MaestroListFile
from matplotlib import pyplot as plt
from pathlib import Path
import numpy as np
from JSB_tools import mpl_hist
from uncertainties import ufloat
from uncertainties import unumpy as unp

#  ====================================================
gamma_Erg = 218
window = 2
max_words = None
t_max = 350
n_overlays = 5
shots = [5, 31, 38, 39, 40]
plot_spectrum = True
day = 'tuesday'
bg_window_offset = 0.5
bin_reduce = 2
sub_bg = False
#  ====================================================
_path = Path.cwd()/'exp_data'
n_pulses = {}


def est_half_live(a):
    a = a[np.argmax(a):]
    return 1.0/np.mean(unp.uarray(a, np.sqrt(a)))


with open(_path/day/'n_pulses') as f:
    lines = f.readlines()
    for shot_num, l in enumerate(lines):
        shot_num += 1
        try:
            n_pulses[shot_num] = int(l)
        except ValueError:
            continue


def get_maestro(num):
    path = _path/day/f'shot{num}.Lis'
    return MaestroListFile(path, max_words=max_words)


def get_cut(a, emin, emax):
    return np.where((a <= emax) & (a >= emin))


ax = None

if plot_spectrum:
    l = get_maestro(13)
    values, bins = np.histogram(l.energies, l.erg_bins)
    mpl_hist(bins, values, np.sqrt(values))

loop_index = 0
for shot_num in shots:
    try:
        pulses = n_pulses[shot_num]
    except KeyError:
        warnings.warn(f"No number of pulses for {shot_num}")
        pulses = np.median(list(n_pulses.values()))
    try:
        l = get_maestro(shot_num)
        same_plot_index = loop_index % n_overlays
        if same_plot_index == 0:
            _, ax = plt.subplots()
            ax.text(.45, 0.9, f'tot. events/t_1/2', transform=ax.transAxes)

        zero_time = np.median(l.realtimes[np.where(l.sample_ready_state != 0)])
        print(f'shot {shot_num} zero time: {zero_time}')

        cut_sig = get_cut(l.energies, gamma_Erg-window//2, gamma_Erg+window//2)
        ergs = l.energies[cut_sig]
        times = l.times - zero_time
        times_sig = l.times[cut_sig]

        values, bins = np.histogram(times_sig, bins=np.linspace(0, t_max, 100))
        bg_values1, _ = np.histogram(times[get_cut(l.energies, gamma_Erg+bg_window_offset+window//2, gamma_Erg+bg_window_offset+window)], bins=np.linspace(0, t_max, 100))
        bg_values2, _ = np.histogram(times[get_cut(l.energies, gamma_Erg-bg_window_offset-window, gamma_Erg-bg_window_offset-window//2)], bins=np.linspace(0, t_max, 100))
        bg_values = (bg_values1 + bg_values2)/2
        values = np.array(values, dtype=np.float)
        if sub_bg:
            values -= bg_values
        values /= pulses
        values *= 600

        # if bin_reduce > 1:
        #     values = np.mean(np.array_split(values, len(values)//bin_reduce), axis=1)

        tot_events = ufloat(sum(values), np.sqrt(sum(values)))

        h = mpl_hist(bins, values, np.sqrt(values), label=f'shot{shot_num}', ax=ax)
        hl = np.log(2)/est_half_live(values)

        ax.legend()
        ax.text(0.45, 0.9-0.1*(same_plot_index+1), f'Shot {shot_num} {tot_events} {hl:.1e}', transform=ax.transAxes)

        loop_index += 1
        # plt.plot(ergs, times_sig)
    except AssertionError:
        print(f"Shot {shot_num} not found")
# l = MaestroListFile('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/tuesday/shot1.Lis')
# l.plot_sample_ready()
# l.plot_count_rate()
plt.show()