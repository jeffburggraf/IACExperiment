import matplotlib.pyplot as plt
import numpy as np
from analysis import Shot, get_merged_time_dependence
# import numpy as np
from JSB_tools import mpl_hist, convolve_gauss
# import matplotlib.pyplot as plt
from JSB_tools.regression import LinearFit
from uncertainties import ufloat
from typing import Dict, List
import uncertainties.unumpy as unp
from JSB_tools import Nuclide, FissionYields, convolve_gauss, fill_between
from lmfit.models import LinearModel


flow = 1.0
first_day = True
b_width =4

shots_warm = {'long': {0.5: [93, 94], 1.0: [96, 97]}, 'short': {0.5: [124, 125], 1.0: [126, 127]}}
shots_cold = {'long': {0.5: [98, 99], 1.0: [100, 101]}, 'short': {0.5: [129, 130], 1.0: [131, 132]}}

time_bins = np.arange(0, 120, 1)

fig, ax = plt.subplots()
line_Styles = {'long': '-', 'short': '--'}
erg = 218.5
for length in ['long', 'short']:
    warm_ys = []
    cold_ys = []
    ls = line_Styles[length]
    for shot_cold, shot_warm in zip(shots_cold[length][flow], shots_warm[length][flow]):
        shot_cold = Shot(shot_cold)
        shot_warm = Shot(shot_warm)
        cold_ys.append(shot_cold.list.get_time_dependence(erg, time_bins, nominal_values=False)[0])
        warm_ys.append(shot_warm.list.get_time_dependence(erg, time_bins, nominal_values=False)[0])
    warm_ys = np.mean(warm_ys, axis=0)
    cold_ys = np.mean(cold_ys, axis=0)

    warm_ys = convolve_gauss(warm_ys, b_width)
    cold_ys = convolve_gauss(cold_ys, b_width)

    warm_ys /= max(warm_ys).n
    cold_ys /= max(cold_ys).n

    # mpl_hist(time_bins, warm_ys, c='red', ax=ax, ls=ls, label=f'Warm; {length}')
    # mpl_hist(time_bins, cold_ys, c='blue', ax=ax, ls=ls, label=f'Cold; {length}')
    fill_between(time_bins, warm_ys, binsxQ=True, ax=ax, ls=ls, color='red', label=f'Warm; {length}')
    fill_between(time_bins, cold_ys, binsxQ=True, ax=ax, ls=ls, c='blue',  label=f'Cold; {length}')
    fig.suptitle("7 cm vs 3 cm filter length")
    ax.set_ylabel("Normalized rate")
    ax.set_xlabel("Time [s]")

plt.show()



