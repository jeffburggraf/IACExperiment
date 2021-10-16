"""




"""
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
from JSB_tools.list_reader import MaestroListFile


# l = MaestroListFile(Path(__file__).parent.parent/'exp_data/Nickel/Nickel.Lis')
# l.plotly(time_bin_width=15*60, time_step=15*60)
# l.pickle()


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
co61 = Nuclide.from_symbol('Co61')
ni57 = Nuclide.from_symbol("Ni57")
u238 = Nuclide.from_symbol('U238')
fiss_xs = u238.gamma_induced_fiss_xs

nickel_spe_chamber = SPEFile(Path.cwd().parent/'exp_data/Nickel/Nickel.Spe')
nickel_spe_vcd: SPEFile = MPA(Path.cwd().parent/'exp_data/Nickel/Nickel.mpa')


def get_obs_nickel57(gamma_line, spe_file, window, efficiency):
    # gamma_line = ni57.decay_gamma_lines[gamma_line_index]
    gamma_i = gamma_line.intensity
    g_erg = gamma_line.erg.n
    # n_gamma_counts = sum(spe_file.get_counts(g_erg-window/2, g_erg+window/2, remove_baseline=False))
    n_gamma_counts, bins = spe_file.get_counts(g_erg-window/2, g_erg+window/2, remove_baseline=True,
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


#  Input efficiency here
ni57_gammaline = ni57.decay_gamma_lines[1]
co61_gammaline = co61.decay_gamma_lines[0]
true_ni_chamber = get_obs_nickel57(ni57_gammaline, nickel_spe_chamber, 4, 0.2625*ufloat(1, 0.0032))
true_co_chamber = get_obs_nickel57(co61_gammaline, nickel_spe_chamber, 4, 0.05*ufloat(1, 0.0032))
#  0.028*ufloat(1, 0.15)
# true_ni_vcd = get_obs_nickel57(1, nickel_spe_vcd, 4, 0.163*ufloat(1, 0.004))
true_ni_vcd = get_obs_nickel57(ni57_gammaline, nickel_spe_vcd, 8, 0.17*ufloat(1, 0.005))


nickel_spe_chamber.plot_erg_spectrum()
ni58_xs = get_ni57_xs(True)
ni62_xs = get_ni57_xs(False)
ni58_xs.plot()
ni62_xs.plot()

out = OutP(Path.cwd()/'sims'/'0_inp'/'outp')

out.get_tally('VCD nickel')
tally_chamber = out.get_tally('Chamber target')  # out.get_cell_by_name('Chamber target').get_tally()
vcd_ni_tally = out.get_tally('VCD nickel')
vcd_cell_tally = out.get_tally('VCD cell')
src_verify_tally = out.get_tally('Source tally')


# calc # of Ni produced according to MCNP
ni_per_src_chamber = np.sum(tally_chamber.dx_per_src * ni58_xs.interp(tally_chamber.energies)) * \
                     (0.68*tally_chamber.cell.atom_density)
co61_per_src_chamber = np.sum(tally_chamber.dx_per_src * ni62_xs.interp(tally_chamber.energies)) * \
                     (0.68*tally_chamber.cell.atom_density)
ni_per_src_vcd = np.sum(vcd_ni_tally.dx_per_src * ni58_xs.interp(vcd_ni_tally.energies)) * \
                 (0.68*vcd_ni_tally.cell.atom_density)

tot_ni_chamber_mcnp = n_electrons * ni_per_src_chamber
tot_co_chamber_mcnp = n_electrons * co61_per_src_chamber
tot_ni_vcd_mcnp = n_electrons * ni_per_src_vcd

#  3/20 to adjust for the fact that there are usually 3 seconds of pulses (n_electrons is assuming 20 seconds)
tally_chamber.plot(norm=n_electrons*3/20*np.pi*chamber_target.radius**2, track_length=False)
vcd_ni_tally.plot(norm=n_electrons*3/20*(vcd_nickel.cross_section_area('z')), track_length=False)

if tot_ni_vcd_mcnp != 0:
    ax = vcd_cell_tally.plot(norm=n_electrons * 3 / 20 * true_ni_chamber / tot_ni_chamber_mcnp * ufloat(1, 0.0), track_length=False,
                             ylabel='Gamma flux [1/cm^2]')

src_verify_tally.plot(track_length=False)


print("Erg [MeV] ; Flux ; flux_err")
for e, f in zip(vcd_cell_tally.energies, vcd_cell_tally.fluxes):
    try:
        f *= n_electrons*3/20*true_ni_vcd/tot_ni_vcd_mcnp*ufloat(1, 0.35)
    except ZeroDivisionError:
        pass
    print(f"{e:.2e}, {f.n:.2e}, {f.std_dev:.2e}")
# ax.sey_ylabel("Particles/")


print(f'n electrons: {n_electrons:.2e}')


s = f"""
Number of Ni57 activation:
         ratio:        {tot_ni_chamber_mcnp/true_ni_chamber}                  {tot_ni_vcd_mcnp/true_ni_vcd}
                      Chamber                  Downstream 
         Model    {tot_ni_chamber_mcnp:.3e}     {tot_ni_vcd_mcnp:.3e}  {tot_ni_chamber_mcnp/tot_ni_vcd_mcnp}
    Gamma spec     {true_ni_chamber:.3e}      {true_ni_vcd:.3}

Number of Co61 activation:
         ratio:        {tot_co_chamber_mcnp/true_co_chamber}               
                      Chamber                 
         Model    {tot_co_chamber_mcnp:.3e}    
    Gamma spec     {true_co_chamber:.3e}     
"""
print(s)

# print('Model # Ni in chamber:', tot_ni_chamber)
# print(' Obs. # Ni in chamber:', true_ni_chamber)
#
#
# print('Model # Ni downstream:', tot_ni_vcd)
# print(' Obs. # Ni downstream:', true_ni_vcd)
plt.show()