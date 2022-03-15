from JSB_tools.maestro_reader import MaestroListFile
from analysis import Shot
from JSB_tools import mpl_hist
import numpy as np
from JSB_tools.spe_reader import save_spe
from matplotlib import pyplot as plt
from pathlib import Path


class Fig:
    def __init__(self):
        self.fig, self.ax = None, None
        self.i = 0

    def __iter__(self):
        return self.gen()

    def gen(self):
        while True:
            if self.fig is None or self.i == len(self.ax):
                self.fig, self.ax = plt.subplots(2, 4, figsize=(10, 10), sharey='all', sharex='all')
                plt.subplots_adjust(hspace=0)
                self.ax = self.ax.flatten()
                self.i = 0

            yield from [self.ax[self.i] for _ in range(3)]

            self.i += 1


#  ====================================================
erg_spectra_min = 60
erg_spectra_max = 2200
erg_spectra_time_max = 300
shot_groups = [list(range(1, 44)),     # 0
               list(range(44, 83)),    # 1
               list(range(107, 119)),  # 2
               list(range(119, 129)),  # 3
               list(range(129, 132)),  # 4
              ]
individual = list(range(98, 107)) + [132, 133, 134, 135, 137, 138, 139]
# lines = [93.63, 137.7, 174.97, 218.59, 511, 602.36, 697.05, 973.9, 1427.7, 1460.8]
lines = [93.63,
         174.97,
         218.59,
         306.83,
         397.44,
         590.24,
         602.36,
         697.05,
         973.9,
         1427.7,
         1460.8,
         1564.64]

#  if cold: 455.49 - Xe-137
#           1118.73  - Kr90
#           121.79  -  Kr90
#           831.69 and 1032.00 -  Rb90
i = 0
plotly = False
time_dep_erg = None
time_dep_sigma = 3
time_dep_b_width = 5
time_step = 5
time_bin_width = 10
plot=False
group_index = 0
#  ====================================================


cwd = Path(__file__).parent

group = list(filter(lambda x: x not in Shot.bad_shots, shot_groups[group_index]))
fig = Fig()
print(group)

l = Shot.sum_shots_list(group[:2])
print( cwd/'fake_spes'/f"group {group_index + 1}")
spe = l.SPE
spe.description = 'Shots ' + str(group)
save_spe(l.SPE, cwd/'fake_spes'/f"group {group_index + 1}")


for k in individual:
    _shot = Shot(k)
    spe = _shot.list.SPE
    spe.description = f'Shot {k}'
    save_spe(spe, cwd/'fake_spes'/f"Shot {k}")


if plotly:
    def f(b0, b1):
        n = Shot.n_shots_array(group, [b0, b1])[0]
        print(b0, b1, n)

        return 1.0/n
    l.plotly(erg_min=60, erg_max=2200, time_step=time_step, time_bin_width=time_bin_width, time_scale_func=f,
             convolve_overlay_sigma=1, remove_baseline=True)
if plot:
    if time_dep_erg is not None:
        time_bins = np.arange(0, l.times[-1], time_dep_b_width)

        sig, _, _ = l.get_time_dependence(time_dep_erg, time_bins, signal_window_kev=time_dep_sigma)
        sig = sig/Shot.n_shots_array(group, time_bins)
        mpl_hist(time_bins, sig)

    ys = []
    axs = set()
    for shot, ax in zip(group, fig):
        axs.add(ax)
        shot = Shot(shot)
        y, bins = shot.list.get_erg_spectrum(erg_spectra_min, erg_spectra_max, time_max=erg_spectra_time_max,
                                             return_bin_edges=True, nominal_values=True, remove_baseline=True)
        ys.append(y)
        i += 1
        mpl_hist(bins, y, ax=ax, label=f"{shot.shotnum}")
    mean_y = np.mean(ys, axis=0)

    for ax in axs:
        mpl_hist(bins, mean_y, label='G_mean', ax=ax, color='black', alpha=0.5, linewidth=5)
        for line in lines:
            ax.axvline(line)
    plt.show()

