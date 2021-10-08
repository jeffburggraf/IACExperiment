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
from inp import vcd_cell, vcd_nickel, chamber_target


def get_ni57_xs():
    with open(Path.cwd()/'Ni57_xs') as f:
        _ergs = []
        _xss = []
        for line in f.readlines():
            erg, xs = line.split()
            erg = float(erg)
            xs = float(xs)*1E-3
            _ergs.append(erg)
            _xss.append(xs)
    return CrossSection1D(_ergs, _xss, 'Ni58(G,N) -> Ni57', 'gamma')


def get_ni57_xs_2():
    with open(Path.cwd()/'Ni57_xs_2') as f:
        _ergs = []
        _xss = []
        for line in f.readlines():
            erg, xs, err_ = line.split(';')
            erg = float(re.match(' *([0-9.]+)', erg).groups()[0])
            xs = float(re.match(' *([0-9.]+)', xs).groups()[0])*1E-3
            _ergs.append(erg)
            _xss.append(xs)
    return CrossSection1D(_ergs, _xss, 'Ni58(G,N) -> Ni57', 'gamma')

# nc_per_second = 60000
# n_seconds = 3  # 60000 per second.
# #  =======================================================================
# nc_per_pulse = nc_per_second*n_seconds
# c_per_pulse = nc_per_pulse * 1E-9
# charge_per_electron = 1.602E-19
# n_electrons = c_per_pulse / charge_per_electron


c_per_second = (192/3.0)*1E-6
charge_per_electron = 1.602E-19
n_electrons = 20*c_per_second / charge_per_electron

ni58 = Nuclide.from_symbol('Ni58')
ni57 = Nuclide.from_symbol("Ni57")
u238 = Nuclide.from_symbol('U238')
fiss_xs = u238.gamma_induced_fiss_xs

nickel_spe_chamber = SPEFile(Path.cwd().parent/'exp_data/Nickel/Nickel.Spe')
nickel_spe_vcd: SPEFile = MPA(Path.cwd().parent/'exp_data/Nickel/Nickel.mpa')


def get_obs_nickel57(gamma_line_index, spe_file, window, efficiency):
    gamma_line = ni57.decay_gamma_lines[gamma_line_index]
    gamma_i = gamma_line.intensity
    g_erg = gamma_line.erg.n
    # n_gamma_counts = sum(spe_file.get_counts(g_erg-window/2, g_erg+window/2, remove_baseline=False))
    n_gamma_counts, bins = spe_file.get_counts(g_erg-window/2, g_erg+window/2, remove_baseline=False,
                                               return_bin_edges=True, deadtime_corr=True)
    n_gamma_counts -= spe_file.deadtime_corr*spe_file.get_baseline_median(g_erg - window / 2, g_erg + window / 2,
                                                                          window_kev=30)
    mpl_hist(bins, n_gamma_counts, title='Bg removed counts. ')
    n_gamma_counts = np.sum(n_gamma_counts)
    frac_decayed_corr = 1.0/(0.5**(spe_file.realtime/ni57.half_life))
    n_ni57 = frac_decayed_corr * n_gamma_counts / gamma_i / efficiency
    return n_ni57


#  Input efficiency here
true_ni_chamber = get_obs_nickel57(1, nickel_spe_chamber, 4, 0.2734*ufloat(1, 0.0032))
#  0.028*ufloat(1, 0.15)
# true_ni_vcd = get_obs_nickel57(1, nickel_spe_vcd, 4, 0.163*ufloat(1, 0.004))
true_ni_vcd = get_obs_nickel57(1, nickel_spe_vcd, 8, 0.17*ufloat(1, 0.005))


nickel_spe_chamber.plot_erg_spectrum()
ni_xs = get_ni57_xs()
# ni_xs = get_ni57_xs_2()
ni_xs.plot()

out = OutP(Path.cwd()/'sims'/'0_inp'/'outp')

# out.get_tally('Chamber target')
out.get_tally('VCD nickel')
tally_chamber = out.get_tally('Chamber target')  # out.get_cell_by_name('Chamber target').get_tally()
vcd_ni_tally = out.get_tally('VCD nickel')
vcd_cell_tally = out.get_tally('VCD cell')


ni_per_src_chamber = 0.68 * np.sum(tally_chamber.dx_per_src * ni_xs.interp(tally_chamber.energies)) * \
                     tally_chamber.cell.atom_density
ni_per_src_vcd = 0.68 * np.sum(vcd_ni_tally.dx_per_src * ni_xs.interp(vcd_ni_tally.energies)) * \
                 vcd_ni_tally.cell.atom_density
tot_ni_chamber_mcnp = n_electrons * ni_per_src_chamber
tot_ni_vcd_mcnp = n_electrons * ni_per_src_vcd

#  3/20 to adjust for the fact that there are usually 3 seconds of pulses (n_electrons is assuming 20 seconds)
tally_chamber.plot(norm=n_electrons*3/20*np.pi*chamber_target.radius**2, track_length=False)
vcd_ni_tally.plot(norm=n_electrons*3/20*(vcd_nickel.cross_section_area('z')), track_length=False)
ax = vcd_cell_tally.plot(norm=n_electrons*3/20*true_ni_vcd/tot_ni_vcd_mcnp*ufloat(1, 0.0), track_length=False,
                         ylabel='Gamma flux [1/cm^2]')

print("Erg [MeV] ; Flux ; flux_err")
for e, f in zip(vcd_cell_tally.energies, vcd_cell_tally.fluxes):
    f *= n_electrons*3/20*true_ni_vcd/tot_ni_vcd_mcnp*ufloat(1, 0.35)
    print(f"{e:.2e}, {f.n:.2e}, {f.std_dev:.2e}")
# ax.sey_ylabel("Particles/")


print(f'n electrons: {n_electrons:.2e}')


s = f"""
Number of Ni activation:
                      Chamber                  Downstream 
         Model    {tot_ni_chamber_mcnp:.3e}     {tot_ni_vcd_mcnp:.3e}
Gamma analysis     {true_ni_chamber:.3e}      {true_ni_vcd:.3}
"""
print(s)
# print('Model # Ni in chamber:', tot_ni_chamber)
# print(' Obs. # Ni in chamber:', true_ni_chamber)
#
#
# print('Model # Ni downstream:', tot_ni_vcd)
# print(' Obs. # Ni downstream:', true_ni_vcd)
plt.show()