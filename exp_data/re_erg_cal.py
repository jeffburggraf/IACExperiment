from JSB_tools.list_reader import MaestroListFile
from analysis import Shot
from JSB_tools import mpl_hist
import numpy as np
from matplotlib import pyplot as plt


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

            yield from [self.ax[self.i] for _ in range(5)]

            self.i += 1


#  ====================================================
erg_spectra_min = 60
erg_spectra_max = 2200
erg_spectra_time_max = 300
shot_groups = [list(range(1, 45))]
lines = [137.72, 174.97, 218.7, 511, 602.36, 1427.7]
i = 0
plotly = True
time_dep_erg = None
time_dep_sigma = 3
time_dep_b_width = 5
time_step = 5
time_bin_width = 10
#  ====================================================


for group in shot_groups:
    group = list(filter(lambda x: x not in Shot.bad_shots, group))
    fig = Fig()
    print(group)

    l = Shot.sum_shots_list(group)

    if plotly:
        def f(b0, b1):
            n = Shot.n_shots_array(group, [b0, b1])[0]
            print(b0, b1, n)

            return 1.0/n
        l.plotly(erg_min=60, erg_max=2200, time_step=time_step, time_bin_width=time_bin_width, time_scale_func=f,
                 convolve_overlay_sigma=1, remove_baseline=True)

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

