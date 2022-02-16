import numpy as np
from matplotlib import pyplot as plt
from JSB_tools.PHITS_tools.PHITStoROOT import PHITSTTreeHelper
from JSB_tools import FileManager, mpl_hist, ROOT_loop, TBrowser
from pathlib import Path
from JSB_tools.TH1 import TH1F
import re
from uncertainties import ufloat
from uncertainties import unumpy as unp

#  =======================================================
sim_dir = 'varygas0.1'  # None for default
# ========================================================

path = Path(__file__).parent

if sim_dir is not None:
    path /= sim_dir

fman = FileManager(path)


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


if sim_dir is None:
    dedx_scale = 1.0
else:
    dedx_scale = float(re.match("vary([0-9.]+)", sim_dir).groups()[0])

print(f'Generating plots of simulations at path {sim_dir}')


for path, attribs in sorted(fman.find_paths(ptrac=True).items(), key=lambda k_v: k_v[1]['my_thickness']):
    tree_helper = PHITSTTreeHelper(path)
    my_thickness = attribs['my_thickness']
    xs.append(my_thickness)

    n_src = {}  # eg {'Xe139': 1100, ....}
    n_stop = {}  # eg {'Xe139': 110, ....},  but then becomes eg {'Xe139': 0.1, ....}  (fraction that stops)

    for evt in tree_helper:
        if evt.is_src:
            try:
                n_src[evt.nuclide.name] += 1
            except KeyError:
                n_src[evt.nuclide.name] = 1
                n_stop[evt.nuclide.name] = {'upstream': 0, 'downstream': 0}
        if evt.is_term and evt.cell == 12 and evt.is_nucleus:
            z = evt.z
            r = np.sqrt(evt.x**2 + evt.y**2)
            ff = evt.nuclide.name
            try:
                stop_data[ff]['data'][my_thickness][0].append(z)
                stop_data[ff]['data'][my_thickness][1].append(r)
            except KeyError:
                fig, (ax_stop_z, ax_stop_r) = plt.subplots(1, 2)
                fig.suptitle(f'{ff}; dedx_sale = {dedx_scale}')
                stop_data[ff] = {'axs': (ax_stop_z, ax_stop_r), 'data': {k: ([z], [r]) for k in [0, 2.5, 5, 10, 20]}}

            if evt.z < 1:
                n_stop[ff]['upstream'] += 1
            else:
                n_stop[ff]['downstream'] += 1

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


fig, ax = plt.subplots()  # fraction stopped vs My thickness

ax.set_xlabel(r"Mylar thickness [$\mu$m]")
ax.set_ylabel(r"Fraction stopped in gas")

for k in ys_up.keys():
    p = ax.errorbar(xs, unp.nominal_values(ys_down[k]), yerr=unp.std_devs(ys_down[k]),
                    label=f'{k} - down')
    color = p[0].get_color()
    tot = np.array(ys_up[k]) + np.array(ys_down[k])
    ax.errorbar(xs, unp.nominal_values(tot), yerr=unp.std_devs(tot),
                label=f'{k} - total', c=color, ls=':')

ax.legend()


plt.show()
