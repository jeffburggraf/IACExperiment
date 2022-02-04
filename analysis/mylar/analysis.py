import numpy as np
from matplotlib import pyplot as plt
from JSB_tools.PHITS_tools.PHITStoROOT import PHITSTTreeHelper
from JSB_tools import FileManager, mpl_hist, ROOT_loop, TBrowser
from pathlib import Path
from JSB_tools.TH1 import TH1F
from uncertainties import ufloat
from uncertainties import unumpy as unp

fman = FileManager(Path(__file__).parent/'sims')
# TBrowser()


xs = []
ys_down = {}
ys_up = {}


for path, attribs in sorted(fman.find_paths(ptrac=True).items(), key=lambda k_v: k_v[1]['my_thickness']):
    tree_helper = PHITSTTreeHelper(path)
    my_thickness = attribs['my_thickness']
    xs.append(my_thickness)

    # hist_src = TH1F(1-2E-4, 1+14E-4, 10)
    # hist_src.Project(tree_helper.tree, "z", "is_src==1")
    # hist_src.Project(tree_helper.tree, '', 'charge_state==54 && is_term==1 && cell')
    n_src = {}
    n_stop = {}
    for evt in tree_helper:
        if evt.is_src:
            try:
                n_src[evt.nuclide.name] += 1
            except KeyError:
                n_src[evt.nuclide.name] = 1
                n_stop[evt.nuclide.name] = {'upstream': 0, 'downstream': 0}
        if evt.is_term and evt.cell == 12:
            if evt.z < 1:
                n_stop[evt.nuclide.name]['upstream'] += 1
            else:
                n_stop[evt.nuclide.name]['downstream'] += 1

    for k, v in n_src.items():
        n_stop[k]['upstream'] = ufloat(n_stop[k]['upstream'], np.sqrt(n_stop[k]['upstream']))/v
        n_stop[k]['downstream'] = ufloat(n_stop[k]['downstream'], np.sqrt(n_stop[k]['downstream']))/v

    for k, v in n_stop.items():
        down = v['downstream']
        up = v['upstream']
        try:
            ys_down[k].append(down)
        except KeyError:
            ys_down[k] = [down]
            ys_up[k] = []

        ys_up[k].append(up)


fig, ax = plt.subplots()

for k in ys_up.keys():
    p = ax.errorbar(xs, unp.nominal_values(ys_down[k]), yerr=unp.std_devs(ys_down[k]),
                    label=f'{k} - down')
    color = p[0].get_color()
    tot = np.array(ys_up[k]) + np.array(ys_down[k])
    ax.errorbar(xs, unp.nominal_values(tot), yerr=unp.std_devs(tot),
                label=f'{k} - total', c=color, ls=':')

ax.legend()


plt.show()
