import numpy as np
from matplotlib import pyplot as plt
from JSB_tools.PHITS_tools.PHITStoROOT import PHITSTTreeHelper
from JSB_tools import FileManager, mpl_hist, ROOT_loop, TBrowser, mpl_hist_from_data
from pathlib import Path
from JSB_tools.TH1 import TH1F
import re
from uncertainties import ufloat
from uncertainties import unumpy as unp

#  =======================================================
sim_dir = 'varygas0.5my0.5'  # None for default
max_n_events = None  # None for None
# ========================================================

path = Path(__file__).parent

if sim_dir is not None:
    path /= sim_dir

fman = FileManager()

varried_mats = []
dedx_scales = []

if sim_dir is not None:
    for m, s in fman.all_files[Path.cwd()/sim_dir].items():
        varried_mats.append(m)
        dedx_scales.append(s)

    fig_title = ', '.join([f"{m}*={s}" for m, s in zip(varried_mats, dedx_scales)])
else:
    fig_title = 'No perturbation'


# TBrowser()
stop_r_bins = np.linspace(0, 2.5, 6)
stop_z_bins = np.linspace(-1.1, 6.5, 20)

# r_weights = 1.0/np.array([rhigh**2 - rlow**2 for rlow, rhigh in zip(stop_r_bins[:-1], stop_r_bins[1:])])

xs = []
ys_down = {}
ys_up = {}

stop_rs = {}
stop_zs = {}
stop_data = {}

where_stop_dicts = {}  #

if sim_dir is not None:
    fman = FileManager(fman.root_directory/sim_dir)

print(f'Generating plots of simulations at path {sim_dir}')

thicknesses_ = np.array([0, 2.5, 5, 10, 20])
ergs_when_leaving_du = {}
prob_stop_vs_erg = {}


for path, attribs in sorted(fman.find_paths(ptrac=True).items(), key=lambda k_v: k_v[1]['my_thickness']):
    tree_helper = PHITSTTreeHelper(path)
    my_thickness = attribs['my_thickness']
    xs.append(my_thickness)

    prob_stop_vs_erg[my_thickness] = {}

    n_events = 0
    n_src = {}  # eg {'Xe139': 1100, ....}
    n_stop = {}  # eg {'Xe139': 110, ....},  but then becomes eg {'Xe139': 0.1, ....}  (fraction that stops)
    th_index = list(thicknesses_).index(my_thickness)

    has_left_du = False
    exit_erg = None

    for evt in tree_helper:
        if evt.is_boundary_crossing and not has_left_du and evt.is_nucleus:
            has_left_du = True
            ff = evt.nuclide.name
            exit_erg = evt.energy
            try:
                ergs_when_leaving_du[ff].append(exit_erg)
            except KeyError:
                ergs_when_leaving_du[ff] = [exit_erg]

        if evt.is_src:
            has_left_du = False
            exit_erg = None

            try:
                n_src[evt.nuclide.name] += 1
            except KeyError:
                n_src[evt.nuclide.name] = 1
                n_stop[evt.nuclide.name] = {'upstream': 0, 'downstream': 0}

        if evt.is_term and evt.is_nucleus:
            ff = evt.nuclide.name

            if ff not in where_stop_dicts:
                where_stop_dicts[ff] = {'my': [0]*5, 'gas': [0]*5, 'chamber': [0]*5}

            if evt.cell == 11:
                where_stop_dicts[ff]['my'][th_index] += 1

            elif evt.term == 'escape':
                where_stop_dicts[ff]['chamber'][th_index] += 1

            elif evt.cell == 12:
                where_stop_dicts[ff]['gas'][th_index] += 1

                if exit_erg is not None:
                    try:
                        prob_stop_vs_erg[my_thickness][ff].append(exit_erg)
                    except KeyError:
                        prob_stop_vs_erg[my_thickness][ff] = [exit_erg]

                z = evt.z
                r = np.sqrt(evt.x**2 + evt.y**2)
                try:
                    stop_data[ff]['data'][my_thickness][0].append(z)
                    stop_data[ff]['data'][my_thickness][1].append(r)
                except KeyError:
                    fig, (ax_stop_z, ax_stop_r) = plt.subplots(1, 2)
                    fig.suptitle(f'{ff}; ({fig_title})')
                    stop_data[ff] = {'axs': (ax_stop_z, ax_stop_r), 'data': {k: ([z], [r]) for k in [0, 2.5, 5, 10, 20]}}

                if evt.z < 1:
                    n_stop[ff]['upstream'] += 1
                else:
                    n_stop[ff]['downstream'] += 1

        n_events += 1
        if max_n_events is not None and n_events > max_n_events:
            break

    for k, v in n_src.items():
        n_stop[k]['upstream'] = ufloat(n_stop[k]['upstream'], np.sqrt(n_stop[k]['upstream']))/v
        n_stop[k]['downstream'] = ufloat(n_stop[k]['downstream'], np.sqrt(n_stop[k]['downstream']))/v

    for k, v in n_stop.items():
        down = v['downstream']
        up = v['upstream']
        try:
            ys_down[k].append(down)
            ys_up[k].append(up)
        except KeyError:
            ys_down[k] = [down]
            ys_up[k] = [up]


fig, axs = plt.subplots(1, 2)
fig.suptitle("Energy of FFs when FFs exit DU foil")
exit_erg_data = {}

for ax, (k, v) in zip(axs, ergs_when_leaving_du.items()):
    y, _ = mpl_hist_from_data(np.linspace(0, max(v), 20), v, title=k, ax=ax)
    exit_erg_data[k] = y
    ax.set_xlabel("Energy [MeV]")
    ax.set_ylabel("counts")

#
# fig, axs = plt.subplots(5, 2)
# fig.suptitle("DU exit energies for FFs that stop in gas")
# for index, (_axs, (th, d)) in enumerate(zip(axs, prob_stop_vs_erg.items())):
#     for ax, (ff, v) in zip(_axs, d.items()):
#         y, bins = np.histogram(v, bins=np.linspace(0, max(v), 20))
#         y = 5*unp.uarray(y, np.sqrt(y))
#         y /= exit_erg_data[ff]
#         mpl_hist(bins, y, ax=ax, title=ff if index == 0 else '')
#         # mpl_hist_from_data(np.linspace(0, max(v), 20), v, title=ff if index == 0 else '', ax=ax)
#         # ax.set_xlabel("Energy [MeV]")
#         # ax.set_ylabel("counts")

for ff, _dd in stop_data.items():
    ax_stop_z, ax_stop_r = _dd['axs']

    ax_stop_z.set_ylabel('counts')
    ax_stop_r.set_ylabel('counts [/cm$^3$]')

    ax_stop_z.set_xlabel('Axial stopping position [cm]')
    ax_stop_r.set_xlabel('Radial stopping position [cm]')

    for thickness, data in _dd['data'].items():
        stop_zs_y, _ = np.histogram(data[0], bins=stop_z_bins)
        stop_rs_y, _ = np.histogram(data[1], bins=stop_r_bins)

        mpl_hist(stop_z_bins, stop_zs_y, label=f'{thickness}', poisson_errors=True, ax=ax_stop_z)
        mpl_hist(stop_r_bins, stop_rs_y, label=f'{thickness}', poisson_errors=True, ax=ax_stop_r)

    ax_stop_z.axvline(1, label='U foil pos.', c='black')

    ax_stop_z.legend(title='Mylar thickness [um]')
    ax_stop_r.legend(title='Mylar thickness [um]')


bar_w = 0.25
bar_x = np.arange(0, len(thicknesses_))
fig, axs = plt.subplots(1, 2)
fig.suptitle(f'DeDx scaling: ({fig_title})')

lines = None
for ax, (ff, d) in zip(axs, where_stop_dicts.items()):
    ax.set_title(ff)
    mult = [-1, 0, 1]
    if lines is None:
        lines = []
        labels = list(d.keys())
    for s, (label, ys) in zip(mult, d.items()):
        alpha = 1 if label == 'gas' else 0.85
        _ = ax.bar(bar_x + s*bar_w, ys, yerr=np.sqrt(ys), ecolor='black', width=bar_w, label=label, capsize=5,
                   alpha=alpha)
        lines.append(_)

    ax.set_xticks(bar_x, list(map(lambda x: fr"{x} $\mu$m", thicknesses_)))
    ax.set_xlabel("Mylar thickness")
    ax.set_ylabel("Counts")
    # ax.legend(title='Where FFs stopped')

fig.legend(lines, labels, title='Where FFs stopped')


fig, ax = plt.subplots()  # fraction stopped vs My thickness

ax.set_xlabel(r"Mylar thickness [$\mu$m]")
ax.set_ylabel(r"Fraction stopped in gas")
ax.set_title(fig_title)


for k in ys_up.keys():
    p = ax.errorbar(xs, unp.nominal_values(ys_down[k]), yerr=unp.std_devs(ys_down[k]),
                    label=f'{k} - down')
    color = p[0].get_color()
    tot = np.array(ys_up[k]) + np.array(ys_down[k])
    ax.errorbar(xs, unp.nominal_values(tot), yerr=unp.std_devs(tot),
                label=f'{k} - total', c=color, ls=':')

ax.legend()

print("Finished!")
plt.show()
