# import re
# import warnings
# import numpy as np
# from JSB_tools.spe_reader import  SPEFile
from mpant_reader import MPA
import matplotlib.pyplot as plt
from pathlib import Path
import datetime
import numpy as np
from JSB_tools import mpl_hist
from JSB_tools.spe_reader import SPEFile
from cal_sources import CalSource
from JSB_tools.nuke_data_tools.gamma_coince import Levels
from JSB_tools.nuke_data_tools import FissionYields
from analysis import Shot
from JSB_tools.nuke_data_tools.gamma_spec import gamma_search
#
# yields = FissionYields('U238', 'gamma', [5], )
#
#
# def weight(n):
#     return np.mean(yields.get_yield(n))
#
#
# for g in gamma_search(296.5, 2, 15, 400, nuclide_weighting_function=weight):
#     print(g)


def _plot(l, plot_time=True, normalize=False):
    fig_erg, ax_erg = plt.subplots()
    if plot_time:
        fig_time, ax_time = plt.subplots()
    time_bins = np.arange(0, 250, 4)

    kwargs_ergs = {'time_max': 310, 'time_min': 0, 'remove_baseline': True, 'erg_min': 200, 'erg_max': 250}
    kwargs_time = {'energy': 218.5, 'bins': time_bins}

    for shot_num in l:
        c = "red" if shot_num != 134 else 'blue'
        shot = Shot(shot_num)
        label = shot.__repr__(['he_flow', 'ar_flow', 'shotnum'])
        y_erg, erg_bins = shot.list.get_erg_spectrum(**kwargs_ergs, return_bins=True)
        if normalize:
            y_erg /= np.max(y_erg)
        mpl_hist(erg_bins, y_erg, label=label, ax=ax_erg, c=c)

        if plot_time:
            rates, _, _ = shot.list.get_time_dependence(**kwargs_time)
            if normalize:
                rates /= np.sum(rates)
            mpl_hist(time_bins, rates, label=label, ax=ax_time)
    if plot_time:
        ax_time.set_xlabel("Time [s]")
        ax_time.set_ylabel("Count rate [arb. units]")

    ax_erg.set_xlabel("Energy")
    ax_erg.set_ylabel("Normalized counts")


_plot(list([120, 133, 129, 131, 134]),  plot_time=False, normalize=True)
_plot(list([133, 129, 131, 134]),  plot_time=True, normalize=True)
# _plot([98, 99, 129, 130])
# Shot(i).list.plotly(remove_background=True)

plt.show()