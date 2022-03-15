import numpy as np
from matplotlib import pyplot as plt
from JSB_tools import mpl_hist, Nuclide
from pathlib import Path
from analysis import Shot
from uncertainties import unumpy as unp

# for s in Shot.find_shots(mylar=20):
#     print(s)


shots = {'100001': {0.25: [[40], [105], [107], [113, 114], [116]], 0.5: [[38], [104], [106], [112], [115]]},
         '010010': {0.25: [[122], [103], [108], [111], [118]], 0.5: [[121], [102], [109], [110], [117]]}}

#  ==============================================================================
ff = 'Sr94'   # Sr94   Sb132'
gamma_index = 0
plotly = True
sigma_erg = 3
#  ==============================================================================

n = Nuclide.from_symbol(ff)
gamma_line = n.decay_gamma_lines[gamma_index]

xs = np.array([0, 2.5, 5, 10, 20])  # my thicknesses

l = None


for pattern, d in shots.items():
    for flow, dd in d.items():
        fig, ax = plt.subplots()
        title = f'{ff}; flow pat = {pattern}; flow rate = {flow} L/min'
        fig.suptitle(title)
        ys = []
        loop_shot_list = []
        _max_time = None
        for shot_list in dd:
            _ = []
            for shot_num in shot_list:
                shot = Shot(shot_num)
                _.append(shot)
                t = shot.list.times[-1]
                if _max_time is None:
                    _max_time = t
                else:
                    if t < _max_time:
                        _max_time = t
            loop_shot_list.append(_)
        print(f"Max time used for {title}: {_max_time}" )
        for index, shot_list in enumerate(loop_shot_list):
            x = xs[index]
            shot_nums_used = []
            list_file = None
            for shot in shot_list:
                print(f"Max time for shot {shot.shotnum}: {shot.list.times[-1]}")
                shot_nums_used.append(shot.shotnum)
                if list_file is None:
                    list_file = shot.list
                else:
                    list_file += shot.list
            counts, bins = list_file.get_erg_spectrum(remove_baseline=True, erg_min=gamma_line.erg.n - sigma_erg,
                                                      erg_max=gamma_line.erg.n + sigma_erg, return_bin_edges=True,
                                                      time_max=_max_time,
                                                      eff_corr=True)
            counts /= len(shot_list)
            ys.append(sum(counts))
            mpl_hist(bins, counts, ax=ax, label=f'{x} um; shots {shot_nums_used}')

        fig = plt.figure()
        fig.suptitle(title)
        plt.errorbar(xs, unp.nominal_values(ys), yerr=unp.std_devs(ys))
        plt.ylim(0)

        plt.show()




# list_file.plotly(convolve_overlay_sigma=1, remove_baseline=True)



