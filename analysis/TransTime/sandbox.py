import numpy as np
from matplotlib import pyplot as plt
from JSB_tools import mpl_hist, TabPlot
from pathlib import Path
from analysis import Shot, MaestroListFile
from JSB_tools import Nuclide
from typing import Tuple
from uncertainties import unumpy as unp

# todo: find nuclei that are not complicated up by parent feeding
# =====================================
nuclide = 'e+'  # use e+ for 511 peak
gamma_indices = [0]

t_max = None
n_time_bins = 20

flow_pat = "100001"

max_shots = 2  # Max number of shots to use for each pipe length

he_flow = 0.25  # also will be Ar flow in L/min

erg_window = 4  # window width in keV for calculating time dependence

time_dep_debug = True  # plot debug in list.get_time_dependence
# =====================================

if nuclide == 'e+':
    nuclide = None
else:
    nuclide = Nuclide.from_symbol(nuclide)

if t_max is None:
    if nuclide is not None:
        t_max = 4*nuclide.half_life.n
    else:
        t_max = 300
    # t_max = min([300, 4*nuclide.half_life.n])

time_bins = np.linspace(0, t_max, n_time_bins)

if nuclide is not None:
    decay_corrections = 0.5**(time_bins[:- 1]/nuclide.half_life.n) - 0.5**(time_bins[1:]/nuclide.half_life.n)
    decay_corrections = 1.0/decay_corrections
else:
    decay_corrections = np.ones(len(time_bins) - 1)

if nuclide is not None:
    print(nuclide, f'half life = {nuclide.half_life.n} seconds')

if nuclide is not None:
    gamma_ergs = [g.erg.n for g in np.array(nuclide.decay_gamma_lines)[gamma_indices]]
else:
    gamma_ergs = [511]

time_dep_tab_plot = TabPlot()


for tube_length in [4.16, 6.14, 9.44, 12.64]:
    print(tube_length)
    shots = []
    for shot in Shot.find_shots(flow=flow_pat, eval_func="'cold' not in self.comment.lower()", flow_stop=0, num_filters=2,
                                llnl_filter_pos=1, tube_len=tube_length, cold_filter=False, mylar=0, beam_duration=3,
                                foil_pos='upstream', he_flow=he_flow, ar_flow=he_flow):
        shots.append(shot)

    if len(shots) > max_shots:
        shots = np.array(shots)[np.random.choice(len(shots), max_shots, False)]

    t_centers = 0.5*(time_bins[1:] + time_bins[:-1])
    rates = None
    y = None
    for shot in shots:
        print(f'\t{shot}')
        gamma_erg = gamma_ergs[0]
        y_i, _, _ = shot.list.get_time_dependence(gamma_erg, time_bins, erg_window, debug_plot=time_dep_debug,
                                                  nominal_values=False, debug_title=str(tube_length))
        if y is None:
            y = y_i
        else:
            y += y_i

    if y is None:
        continue

    y *= decay_corrections

    t_dep_ax = time_dep_tab_plot.new_ax(tube_length)

    mpl_hist(time_bins, y, ax=t_dep_ax)

plt.show()





