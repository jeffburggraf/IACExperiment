import turtle
import random

import matplotlib.pyplot as plt

from analysis import Shot

l_tot = None

fig, axs = plt.subplots(2, 2, sharex='all')

axs = axs.flatten()
ax = axs[0]

sss = range(1, 135, 10)


for index, s in enumerate(sss):
    shot = Shot(s)
    shot.list.energy_binned_times
    shot.list.energy_binned_times
    shot.list.energy_binned_times
    if index %4 == 0:
        ax = axs[index//4]
    if l_tot is None:
        l_tot = shot.list
    else:
        l_tot += shot.list
    shot.list.plot_erg_spectrum(ax=ax, label=shot.shotnum, erg_max=1500)

for ax in axs:
    l_tot.plot_erg_spectrum(ax=ax, label='tot', erg_max=1500, scale=1.0/len(sss))
    ax.legend()

plt.show()

