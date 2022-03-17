import numpy as np
from matplotlib import pyplot as plt
from JSB_tools import mpl_hist
from pathlib import Path
from analysis import Shot, MaestroListFile
from JSB_tools import Nuclide
from typing import Tuple
from uncertainties import unumpy as unp

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 13})

#  =====================================================================================================================
erg_ranges = [(64, 250), (255, 511), (511, 1000), (1026, 1446)]
times_ranges = [(0, 15), (30, 300)]
make_rate=True
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

override_labels = {('Kr90', 121.79): "$^{90}$Kr/$^{99}$Y"}

nuclides_gammalines = {'Xe139': [0, 3, 4, 5, 6, 7],  # See special cases for gamma 1
                       'Nb101': [2, 3],
                       'Rb94': [0],
                       'Nb100': [0, 2],
                       'Xe140': [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 14],
                       "Xe137": [0],
                       'I134': [0],
                       'Kr91': [0, 1],
                       'Zr100': [0, 1],
                       'Nb99': [0, 1],
                       'Kr90': [0, 1, 2],
                       'Kr92': [0, 1],
                       'Kr89': [0, 1],
                       'Sr95': [0],
                       'Sb132': [0, 1, 3],
                       'Xe138': [0, 1],
                       'Cs140': [0],
                       'Sr94': [0],
                       'Rb89': [0, 1],
                       'Rb91': [0],
                       'Zr99': [0, 1, 2],
                       'Rb90': [0],
                       'Ce145': [0],
                       'Mo104': [0],
                       'Nb102': [0],
                       'Ba143': [0],
                       '?': [515.2]
                       }
"""
Verified:
Ba143 @ 210. hl seems a bit longer than 15 s... but what else would it be? 
Kr92 @ 142.3 and 1218 keV (shortest half life observed at 1.8 s)
Kr90/Y99 at 121.79/121.7 keV in the cold and in warm spectra, respectively
Kr91 @ 108 & 506
Nb99 @ 137.7
Rb94 (hl=2.7)
Nb100 (hl=2.3 s) is being fed by parents, giving a longer apparent half life. . 
Nb101 @ 276 (hl=8)
Nb102 (hl=4 s) contributes to the line at 296 (along with Xe139). You can see a shift from 296 to 296.6. 
    You can also track the rel frac of the Xe139 lines. At early times, the 296 line is 2x more than should be rel 
    to the 218 line.

Notes:
Add special case for Xe139 and Ba143. Blue line is Ba143, right? 
Nb100 is shorted, but it is bottle necked by feeding from Zr100
very rapidly falling line at 306.7 might be Cs143, but its not quite the right energy (~.3 keV diff). 
Peak around 258.4 is Xe138 (mostly), Pr149, Sb130, and La146_m1
I137 and Kr92 have close lines near 1218. 
Lines at 724 are a mixture of Zr95 build up and Ce145

Ba141: Expected gamma lines at 277, 304.2, & 190.3 are NOT seen. Models suggest this should be highest yield for t<100 s
"""
#  =====================================================================================================================
# bg = Shot.background_spe()

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
        warm_list.plotly(erg_max=2000, erg_min=50, time_max=150, time_bin_width=5, time_step=1.5, leg_label='Warm',
                         convolve_overlay_sigma=1,
                         remove_baseline=True)
        cold_list.plotly(erg_max=2000, erg_min=50, time_bin_width=7.5, time_step=3, time_max=300, leg_label="Cold",
                         convolve_overlay_sigma=1, remove_baseline=True)

    print(f"{warm_list.n} warm shots used")
    print(f"{cold_list.n} cold shots used")

else:
    warm_list = None
    cold_list = None

fig, axs = plt.subplots(len(erg_ranges), 1, figsize=(14, 10))

plt.subplots_adjust(left=0.1, right=0.9, hspace=0.4)


nuclides_ids = {}

ys_dict = [{'hot': None, 'cold': None} for i in range(len(axs))]  # dict of y values for hot cold and each axis
xs_dict = [{'hot': None, 'cold': None} for i in range(len(axs))]  # dict of x bins for hot cold and each axis


def annotate(_ax, nuclide_name, erg, y_value):
    rotation = 76
    if nuclide_name == '?':
        n_id = "?"
        rotation = 0
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

    _ax.plot([erg, erg], [y_value*1.05, ymax], color='black', linewidth=1)

    _ax.set_ylim(None, ymax)
    a = 0
    if (nuclide_name, erg) in strattler:
        a = strattler[(nuclide_name, erg)]

    try:
        label = override_labels[(nuclide_name, erg)]
    except KeyError:
        if nuclide_name != '?':
            n = Nuclide.from_symbol(nuclide_name)
            label = f"$^{{{n.A}}}${n.atomic_symbol}"
        else:
            label = '?'

    _ax.text(erg + a, ymin + 1.03*(ymax-ymin), label, rotation=rotation, size='medium')


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

strattler = {}  # for strattling pairs of closeby lines.

handles, labels = None, []

for index, erg_range in enumerate(erg_ranges):
    e1, e2 = erg_range
    ax = axs[index]
    max_ys = None

    interp_may_ys_x = None
    is_cold=None

    _handles = []

    for tmin, tmax in times_ranges:
        spec: MaestroListFile
        if tmax < warm_shots_thresh:
            spec = warm_list
            is_cold = False
        else:
            spec = cold_list
            is_cold = True

        if spec is not None:
            y, bins = spec.get_erg_spectrum(e1, e2, time_min=tmin, time_max=tmax, return_bin_edges=True,
                                            remove_baseline=True, nominal_values=True)
            scale = spec.n
        else:
            bins = np.linspace(e1, e2, 101)
            y = tmax*(np.sin(2*np.pi*np.linspace(0, e2-e1, 100)/(e2-e1)) + 1)

            scale = 1

        y /= scale
        if make_rate:
            y /= tmax
        _, lines = mpl_hist(bins, y, ax=ax, return_line=True)

        if handles is None:
            _handles.append(lines[0])
            labels.append(f"${tmin}<t<{tmax}$ seconds")

        ax.set_ylabel("counts/s")
        ax.margins(x=0)
        # ax.ticklabel_format(axis='y')

        if interp_may_ys_x is None:
            interp_may_ys_x = 0.5 * (bins[1:] + bins[:-1])
            max_ys = y
        else:
            new_interp = 0.5 * (bins[1:] + bins[:-1])
            y = np.interp(interp_may_ys_x, new_interp, y)
            max_ys = np.max([max_ys, y], axis=0)

        for i in range(0, len(max_ys), 3):
            i1, i2 = i, i + 3
            i2 = min([len(max_ys) - 1, i2])
            if i1 == i2:
                continue

            max_ys[i1: i2] = max(max_ys[i1: i2])

        if is_cold:
            ys_dict[index]['cold'] = y[:]
            xs_dict[index]['cold'] = bins
        else:
            ys_dict[index]['hot'] = y[:]
            xs_dict[index]['hot'] = bins

    if handles is None:
        handles = _handles

    for i1 in range(len(erg_ranges_dict[erg_range])-1):
        for i2 in range(i1 + 1, len(erg_ranges_dict[erg_range])):
            erg1 = erg_ranges_dict[erg_range][i1][1]
            erg2 = erg_ranges_dict[erg_range][i2][1]

            if abs(erg1 - erg2) < 2:
                print(erg_ranges_dict[erg_range][i1], erg_ranges_dict[erg_range][i2])
                if abs(erg1 - erg2) < 1:
                    dist = 2.5
                else:
                    dist = 1
                strattler[tuple(erg_ranges_dict[erg_range][i1])] = -dist
                strattler[tuple(erg_ranges_dict[erg_range][i2])] = dist

    for nuclide_name, erg in erg_ranges_dict[erg_range]:
        i = np.searchsorted(interp_may_ys_x, erg, side='left')
        annotate_y = max_ys[i]
        annotate(ax, nuclide_name, erg, annotate_y)

    ax.minorticks_on()

axs[-1].set_xlabel("$\gamma$ energy [keV]")
print(handles, labels)
fig.legend(handles, labels, loc='upper right')

# ================= Special Cases ======================
axs_i = 1
h_or_c = 'cold'
n = Nuclide.from_symbol("Xe139")
gamma_index = 1

erg = n.decay_gamma_lines[gamma_index].erg.n
i = np.searchsorted(xs_dict[axs_i][h_or_c], erg, side='right') - 1
y_annotate = 1.06*ys_dict[axs_i][h_or_c][i]
ylim0 = axs[axs_i].get_ylim()
axs[axs_i].text(erg+3, ylim0[1], f"$^{{{n.A}}}${n.atomic_symbol}", rotation=76)
axs[axs_i].plot([erg, erg + 2.5], [y_annotate, ylim0[1]], c='black')
axs[axs_i].set_ylim(ylim0)
# ======================================================


for n_name, indices in nuclides_gammalines.items():
    if n_name == '?':
        continue
    n = Nuclide.from_symbol(n_name)
    print(f"{n_name}; {n.human_friendly_half_life()}")
    gamma_lines = n.decay_gamma_lines
    gamma_lines = [gamma_lines[i] for i in indices]
    for g in gamma_lines:
        print(f"\t{g}")

plt.subplots_adjust(hspace=0.48)
plt.show()


