import re
from JSB_tools.nuke_data_tools import Nuclide
from matplotlib import pyplot as plt
from pathlib import Path
from JSB_tools.MCNP_helper.outp_reader import OutP
import numpy as np
from JSB_tools import mpl_hist
from JSB_tools.spe_reader import SPEFile
from mpant_reader import MPA
from JSB_tools.nuke_data_tools import CrossSection1D
import uncertainties.unumpy as unp
from uncertainties import ufloat
from analysis import Shot
from inp import vcd_cell, vcd_nickel, chamber_target
from JSB_tools.maestro_reader import MaestroListFile


c_per_second = (192/3.0)*1E-6
charge_per_electron = 1.602E-19
n_electrons = 20*c_per_second / charge_per_electron


nickel_spe_chamber = SPEFile(Path.cwd().parent/'exp_data/Nickel/Nickel.Spe')


def get_ni57_xs(ni58_bool):
    if ni58_bool:
        names = ("Ni57(g,n0).csv", "Ni57(g,n1).csv", "Ni57(g,nc).csv")
    else:
        names = ('Ni62(g,p).csv',)
    outs = []
    for name in names:
        xs = []
        ys = []
        yerrs = []
        with open(Path.cwd()/'xs_data'/name) as f:
            for line in f.readlines():
                try:
                    values = tuple(map(float, line.split(';')))
                except ValueError:
                    continue
                try:
                    x, y, yerr = values
                except ValueError:
                    x, y = values
                    yerr = 0
                x *= 1E-6
                xs.append(x)
                ys.append(y)
                yerrs.append(yerr)
        outs.append((xs, ys))

    x_interp = outs[np.argmax(list(map(lambda x: len(x[0]), outs)))][0]
    out = np.sum([np.interp(x_interp, x, y) for x, y in outs], axis=0)

    return CrossSection1D(x_interp, out, 'Ni58(G,N) -> Ni57', 'gamma')

ni57 = Nuclide.from_symbol("Ni57")


def get_obs_nickel57(gamma_line, spe_file, window, efficiency):
    # gamma_line = ni57.decay_gamma_lines[gamma_line_index]
    gamma_i = gamma_line.intensity
    g_erg = gamma_line.erg.n
    # n_gamma_counts = sum(spe_file.get_counts(g_erg-window/2, g_erg+window/2, remove_baseline=False))
    n_gamma_counts, bins = spe_file.get_counts(g_erg-window/2, g_erg+window/2, remove_baseline=False,
                                               return_bin_edges=True, deadtime_corr=True, debug_plot=True,
                                               baseline_method='median', baseline_kwargs={'window_kev': 40})
    # n_gamma_counts -= spe_file.deadtime_corr*spe_file.get_baseline_median(g_erg - window / 2, g_erg + window / 2,
    #                                                                       window_kev=30)
    # mpl_hist(bins, n_gamma_counts, title='Bg removed counts. ')
    n_gamma_counts = np.sum(n_gamma_counts)
    frac_decayed_corr = 1.0/(1-0.5**(spe_file.realtime/ni57.half_life))
    n_ni57 = frac_decayed_corr * n_gamma_counts / gamma_i / efficiency
    print('Deadtime corr:', spe_file.deadtime_corr, gamma_i)
    return n_ni57


ni_true = get_obs_nickel57(ni57.decay_gamma_lines[1], nickel_spe_chamber, 4, 0.2625*ufloat(1, 0.0032))
out = OutP(Path.cwd()/'sims'/'simple_ni'/'outp')
tally = out.get_tally('Target')
ni_atom_density = tally.cell.atom_density

ni_xs = get_ni57_xs(True)


model_ni = sum(tally.dx_per_src*ni_xs.interp(tally.energies)*ni_atom_density*n_electrons)
print('model: ', model_ni)
print('meas: ', ni_true)

tally.plot()
plt.show()
