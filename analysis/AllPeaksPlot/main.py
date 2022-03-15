import numpy as np
from matplotlib import pyplot as plt
from JSB_tools import mpl_hist
from pathlib import Path
from analysis import Shot, MaestroListFile
from JSB_tools import Nuclide
from typing import Tuple
from uncertainties import unumpy as unp

#  =====================================================================================================================
erg_ranges = [(64, 250), (250, 508), (515, 1000), (1000, 1500)]
times_ranges = [(0, 15), (30, 120)]

# A max time (from `times` above) below this will use warm shots (for better statistics).
warm_shots_thresh = 20

unpickle = True

warm_shot_limit = 20  # for speed. Too many warm shots. Do 128 for fastest.

n_warm_shots_max = 45
n_cold_shots_max = 25


warm_list_path = Path(__file__).parent/'warm_list.pickle'
cold_list_path = Path(__file__).parent/'cold_list.pickle'

plotly = False

no_draw = False

nuclides_gammalines = {'Xe139': [0, 1, 3, 4, 5, 6, 7],
                       'Kr90': [0, 1, 2],
                       'Nb101': [2, 3],
                       'Rb94': [0],
                       'Kr89': [0, 1],
                       'Nb100': [0, 2],
                       'Xe140': [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 14, 15, 16],
                       "Xe137": [0],
                       'I134': [0],
                       'Kr91': [0, 1],
                       'Zr100': [0, 1],
                       'Nb99': [0, 1],
                       'Kr92': [0, 1],
                       'Sr95': [0],
                       'Sb132': [0, 1],
                       'Xe138': [0],
                       'I137': [0],
                       'Cs140': [0],
                       'Sr94': [0],
                       'Rb89': [0, 1],
                       'Zr99': [0, 1]
                       }
"""
Kr92 @ 142.3 is the shortest half life observed at 1.8 s
Nb100 is shorted, but it is bottle necked by feeding from Zr100
very rapidly falling line at 306.7 might be Cs143, but its not quite the right energy (~.3 keV diff). 
Peak around 258.4 is Xe138 (mostly), Pr149, Sb130, and La146_m1
I137 and Kr92 have close lines near 1218. 
Lines at 724 are a mixtue of Zr95 build up and Ce145
306 keV unknown

"""
#  =====================================================================================================================

if not no_draw:
    if not unpickle:
        warm_shots = Shot.find_shots(tube_len=4.16, cold_filter=False, mylar=0, flow_stop=0, beam_duration=3,
                                     num_filters=2,
                                     eval_func=f'self.he_flow + self.ar_flow >= 1.0 and self.shotnum>{warm_shot_limit}')[:n_warm_shots_max]

        cold_shots = Shot.find_shots(tube_len=4.16, cold_filter=True)[:n_cold_shots_max]

        warm_list = [None]
        cold_list = [None]

        for l, shots in zip([cold_list, warm_list], [cold_shots, warm_shots]):
            for shot in shots:
                print(shot)
                if l[0] is None:

                    l[0] = shot.list
                    l[0].n = 0
                else:
                    l[0].__add__(shot.list, False, False)

                l[0].n += 1

            print("Warm:")
        warm_list = warm_list[0]
        cold_list = cold_list[0]

        warm_list.pickle(warm_list_path, meta_data={'n': warm_list.n})
        cold_list.pickle(cold_list_path, meta_data={'n': cold_list.n})

    else:
        cold_list = MaestroListFile.from_pickle(cold_list_path)
        warm_list = MaestroListFile.from_pickle(warm_list_path)

    if plotly:
        warm_list.plotly(erg_max=2000, erg_min=50, time_max=100, time_bin_width=5, time_step=1.5, leg_label='Warm',
                         convolve_overlay_sigma=1,
                         remove_baseline=True)
        cold_list.plotly(erg_max=2000, erg_min=50, time_bin_width=10, time_step=5, time_max=300, leg_label="Cold",
                         convolve_overlay_sigma=1, remove_baseline=True)

    print(f"{warm_list.n} warm shots used")
    print(f"{cold_list.n} cold shots used")

else:
    warm_list = None
    cold_list = None

fig, axs = plt.subplots(len(erg_ranges), 1, figsize=(14, 10))

plt.subplots_adjust(left=0.1, right=0.9, hspace=0.4)


nuclides_ids = {}


def annotate(_ax, nuclide_name, erg, y_value):
    if nuclide_name == '?':
        n_id = "?"
    else:
        try:
            n_id = nuclides_ids[nuclide_name]
        except KeyError:
            n_id = 1
            while n_id in nuclides_ids.values():
                n_id += 1
            nuclides_ids[nuclide_name] = n_id

    ymax = _ax.get_ylim()[1]
    ymin = _ax.get_ylim()[0]

    _ax.plot([erg, erg], [y_value, ymax], color='black', linewidth=1.5)

    _ax.set_ylim(None, ymax)

    # _ax.text(erg, ymin + 1.03*(ymax-ymin), str(n_id))
    # plt.text()
    n = Nuclide.from_symbol(nuclide_name)
    _ax.text(erg, ymin + 1.03*(ymax-ymin), f"$^{{{n.A}}}${n.atomic_symbol}", rotation=76)


erg_ranges_dict = {e1e2: [] for e1e2 in erg_ranges}  # {(e1, e2): [[Xe139, 219.59], [Kr90, 221.2], ...]}

for n, gamma_indices in nuclides_gammalines.items():
    if n == '?':
        nuclide_name = "?"
        gamma_ergs = gamma_indices
    else:
        nuclide = Nuclide.from_symbol(n)
        nuclide_name = nuclide.name
        gamma_ergs = [nuclide.decay_gamma_lines[i].erg.n for i in gamma_indices]
    for erg in gamma_ergs:
        for e1, e2 in erg_ranges:
            if e1 <= erg <= e2:
                erg_ranges_dict[(e1, e2)].append([nuclide_name, erg])

for k, v in erg_ranges_dict.items():
    erg_ranges_dict[k] = list(sorted(v, key=lambda x: x[1]))

for index, erg_range in enumerate(erg_ranges):
    e1, e2 = erg_range
    ax = axs[index]
    max_ys = None
    interp_may_ys_x = None

    for tmin, tmax in times_ranges:
        spec: MaestroListFile
        if tmax < warm_shots_thresh:
            spec = warm_list
        else:
            spec = cold_list

        if spec is not None:
            y, bins = spec.get_erg_spectrum(e1, e2, time_min=tmin, time_max=tmax, return_bin_edges=True,
                                            remove_baseline=True, nominal_values=True)
            scale = spec.n
        else:
            bins = np.linspace(e1, e2, 101)
            y = tmax*(np.sin(2*np.pi*np.linspace(0, e2-e1, 100)/(e2-e1)) + 1)

            scale = 1

        y /= (tmax*scale)
        mpl_hist(bins, y, ax=ax, label=f"{tmin}<t<{tmax}")
        ax.set_ylabel("counts")
        ax.margins(x=0)
        ax.ticklabel_format(axis='y', scilimits=(0, 0))

        if interp_may_ys_x is None:
            interp_may_ys_x = 0.5 * (bins[1:] + bins[:-1])
            max_ys = y
        else:
            new_interp = 0.5 * (bins[1:] + bins[:-1])
            y = np.interp(interp_may_ys_x, new_interp, y)
            max_ys = np.max([max_ys, y], axis=0)

    for i1 in range(len(erg_ranges_dict[erg_range])-1):
        for i2 in range(i1 + 1,len(erg_ranges_dict[erg_range])):
            erg1 = erg_ranges_dict[erg_range][i1][1]
            erg2 = erg_ranges_dict[erg_range][i2][1]
            if abs(erg1 - erg2) < 2:
                print(erg_ranges_dict[erg_range][i1], erg_ranges_dict[erg_range][i2])

    for nuclide_name, erg in erg_ranges_dict[erg_range]:
        annotate(ax, nuclide_name, erg, max_ys[np.searchsorted(interp_may_ys_x, erg, side='left')])

plt.show()


