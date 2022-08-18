# import lmfit
import pickle

from lmfit import Model
from lmfit.models import LinearModel
import numpy as np
from matplotlib import pyplot as plt
from JSB_tools import mpl_hist, TabPlot
from pathlib import Path
from analysis import Shot, MaestroListFile
from JSB_tools import Nuclide, convolve_gauss
from typing import Tuple
from uncertainties import unumpy as unp
from uncertainties import ufloat
from scipy.signal import argrelextrema




# todo: find nuclei that are not complicated up by parent feeding
# =====================================
nuclide = 'Xe139'  # use e+ for 511 peak
gamma_indices = [0, 1, 2, 3]

t_max = None
n_time_bins = 50

convolve = 4  # convolve in s. None for no convolve

flow_pat = "010010"  # only really works for 010010 due to lack of data

max_shots = None  # Max number of shots to use for each pipe length

he_flow = 0.25  # also will be Ar flow in L/min

erg_window = 4  # window width in keV for calculating time dependence

time_dep_debug = False  # plot debug in list.get_time_dependence
# =====================================

if nuclide == 'e+':
    nuclide = None
else:
    nuclide = Nuclide.from_symbol(nuclide)

if t_max is None:
    if nuclide is not None:
        t_max = min([350, 6*nuclide.half_life.n])
    else:
        t_max = 350

time_bins = np.linspace(0, t_max, n_time_bins)
t_centers = 0.5 * (time_bins[1:] + time_bins[:-1])

with open(Path(__file__).parent.parent/'ff_decay_rates.pickle', 'rb') as f:
    decay_rates_data = pickle.load(f)
    decay_times = pickle.load(f)

decay_rates = None

for k, v in decay_rates_data[nuclide.name].items():
    if decay_rates is None:
        decay_rates = v
    else:
        decay_rates += v

decay_rates = np.interp(t_centers, decay_times, decay_rates)
decay_rates *= 1.0/sum(decay_rates)
decay_corrections = 1.0/decay_rates

mpl_hist(time_bins, decay_rates)

plt.xlabel("Time [s]")
plt.ylabel("Rel. amount decayed")

# cumsum_decay_rates = np.cumsum(decay_rates)
# decay_corrections = cumsum_decay_rates[1:] - cumsum_decay_rates[:-1]


# if nuclide is not None:
#     decay_corrections = 0.5**(time_bins[:- 1]/nuclide.half_life.n) - 0.5**(time_bins[1:]/nuclide.half_life.n)
#     decay_corrections = 1.0/decay_corrections
# else:
#     decay_corrections = np.ones(len(time_bins) - 1)

if nuclide is not None:
    print(nuclide, f'half life = {nuclide.half_life.n} seconds')

if nuclide is not None:
    gamma_ergs = [g.erg.n for g in np.array(nuclide.decay_gamma_lines)[gamma_indices]]
else:
    gamma_ergs = [511]

time_dep_tab_plot = TabPlot()


def truncate(xbins, y):
    x = 0.5*(xbins[1:] + xbins[:-1])
    yerr = unp.std_devs(y)
    y = unp.nominal_values(y)
    max_i = argrelextrema(y, np.greater)[0]
    if not len(max_i):
        max_i = len(y)
    else:
        max_i = max_i[0]

    x = x[:max_i]
    y = y[:max_i]
    yerr = yerr[:max_i]
    xbins = xbins[:max_i + 1]

    # t = sum(y)
    # y /= t
    # yerr /= t

    return x, (xbins[1:] - xbins[:-1])/(2*np.sqrt(3)), y, yerr


tube_lengths = [4.16, 6.14, 9.44, 12.64]
counts_v_times = []
fig, inv_ax = plt.subplots()

for tube_length in tube_lengths:
    print(tube_length)
    shots = []
    for shot in Shot.find_shots(flow=flow_pat, eval_func="'cold' not in self.comment.lower()", flow_stop=0, num_filters=2,
                                llnl_filter_pos=1, tube_len=tube_length, cold_filter=False, mylar=0, beam_duration=3,
                                foil_pos='upstream', he_flow=he_flow, ar_flow=he_flow):
        shots.append(shot)

    if not len(shots):
        continue

    if max_shots is not None and len(shots) > max_shots:
        shots = np.array(shots)[np.random.choice(len(shots), max_shots, False)]

    rates = None

    max_acq_time = None
    for shot in shots:
        t = shot.list.times[-1]
        if max_acq_time is None:
            max_acq_time = t
        else:
            if t < max_acq_time:
                max_acq_time = t

    if t_max > max_acq_time:
        t_max = max_acq_time

    i = np.searchsorted(time_bins, max_acq_time, side='right') - 1
    i = min([len(time_bins) - 1, i])
    trunc_time_bins = time_bins[:i + 1]
    trunc_decay_corrections = decay_corrections[:i]
    trunc_t_centers = t_centers[:i]

    y = None

    for shot in shots:
        print(f'\t{shot}')
        for gamma_erg in gamma_ergs:
            y_i, _, _ = shot.list.get_time_dependence(gamma_erg, trunc_time_bins, erg_window, debug_plot=time_dep_debug,
                                                      nominal_values=False, debug_title=str(tube_length))
            if y is None:
                y = y_i
            else:
                y += y_i

    if y is None:
        continue

    y_corrected = y*trunc_decay_corrections

    if convolve is not None:
        b_w = time_bins[1] - time_bins[0]
        y_corrected = convolve_gauss(y_corrected, convolve/b_w)
        y = convolve_gauss(y, convolve/b_w)

    assert len(trunc_time_bins) - 1 == len(y_corrected)
    assert len(trunc_time_bins) - 1 == len(y)

    # counts.append(trunc_decay_corrections)
    # times.append(trunc_t_centers)

    t_dep_axs = time_dep_tab_plot.new_ax(tube_length, 1, 2, suptitle=f"tube len = {tube_length}")

    mpl_hist(trunc_time_bins, y_corrected, ax=t_dep_axs[0], label='decay corr.')
    mpl_hist(trunc_time_bins, y, ax=t_dep_axs[1], label='raw decay rate')

    inv_x, inv_xerr, inv_y, inv_yerr = truncate(trunc_time_bins, y_corrected)

    inv_ax.errorbar(inv_x, inv_y, xerr=inv_xerr, yerr=inv_yerr, label=f'{tube_length}')

inv_ax.legend()

t_displacements = []


plt.show()

