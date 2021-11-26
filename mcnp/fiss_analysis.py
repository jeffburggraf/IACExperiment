"""
Thoughts Nov 1:
    Is Kr89 really Kr89? Check with better statistics by summiong cold filter shots
    How to do decay correction?
    Make more options to disable certain corrections.

"""
import warnings
from JSB_tools.nuke_data_tools import Nuclide, DecayedFissionYields
from matplotlib import pyplot as plt
from pathlib import Path
from JSB_tools.MCNP_helper.outp_reader import OutP
import numpy as np
from JSB_tools.spe_reader import SPEFile
from JSB_tools.nuke_data_tools import CrossSection1D
import uncertainties.unumpy as unp
from uncertainties import ufloat
from analysis import Shot
from JSB_tools.nuke_data_tools import FissionYields
from stopping_powers import get_fraction
from lmfit.model import Model

c_per_second = (192/3.0)*1E-6
charge_per_electron = 1.602E-19
n_electrons = 3*c_per_second / charge_per_electron

#  ======================================================
model_correction = 1.0/1.5
nuclide = 'Xe139'  # Kr89, Sb132, Sr94
shot_num = 134
suppress_upstream = True  # Do or don't include FF's which escape from upstream of foil
do_decay_corr = False
assume_trans_time = 20
gamma_index = 0

#  ======================================================

shot = Shot(shot_num)
# shot.list.plotly(remove_baseline=True)


ni_meas_scale = n_electrons * model_correction

ff = Nuclide.from_symbol(nuclide)


def decay_model(x, time_shift, scale, lambda_):
    return scale*np.e**(-(x-time_shift)*lambda_)*lambda_


if do_decay_corr:
    max_time = shot.max_time

    if isinstance(assume_trans_time, (float, int)):
        rise_time = assume_trans_time
    else:
        time_dep, _, time_bins = shot.list.get_time_dependence(ff.decay_gamma_lines[0].erg.n, debug_plot='simple',
                                                                signal_window_kev=1.5)
        shot.list.plot_time_dependence(ff.decay_gamma_lines[0].erg.n, signal_window_kev=1.5)
        rise_time = 0.5*(time_bins[1:] + time_bins[:-1])[np.where(time_dep > 0.5*max(time_dep))[0][0]]

        print(f"rise_time: {rise_time}")
    decay_corr = 1.0 / (0.5 ** (rise_time / ff.half_life) - 0.5 ** (max_time / ff.half_life))

else:
    decay_corr = 1

gamma_line = ff.decay_gamma_lines[gamma_index]
fit_ergs = [gamma_line.erg.n]
if nuclide == 'Xe139':
    fit_ergs.append(221)
elif nuclide == 'Kr89':
    fit_ergs.append(218.5)
fit_result = shot.llnl_spe.multi_peak_fit(fit_ergs, fit_window=30, eff_corr=True, debug_plot=True,
                                          baseline_method='median')
print(fit_result.fit_report())


print(f"Fission fragment: {ff}")
print(f"Gamma line used: {gamma_line}")
print(f'decay_corr: {decay_corr}')

amp_error = fit_result.params['_0amplitude'].stderr
if amp_error is None:
    amp_error = fit_result.params['_0amplitude'].value*0.1
    warnings.warn("Fit failed to calculate errors!")

n_ff_meas = ufloat(fit_result.params['_0amplitude'].value, amp_error) / gamma_line.intensity
n_ff_meas *= decay_corr

# shot.llnl_spe.plot_erg_spectrum()

outp = OutP(Path(__file__).parent/'sims'/'1_inp'/'outp')
tally_up = outp.get_tally('Active up')
tally_down = outp.get_tally('Active down')
# tally_down.plot()
# tally_up.plot()


atom_density = tally_up.cell.atom_density

u238 = Nuclide.from_symbol('U238')
# u238.gamma_induced_fiss_xs.plot()
yields = FissionYields('U238', 'gamma', tally_down.energies, independent_bool=True)
# yields.plot(nuclide)

ff_yield = np.average(yields.get_yield(nuclide), weights=tally_down.nominal_fluxes)


print('Atom density: ', atom_density)
print(f'{nuclide} cumulative yield: ', ff_yield)

# plt.plot(tally_down.energies, u238.gamma_induced_fiss_xs.interp(tally_down.energies))

dedx_sub_nuclide = nuclide  # in case an equiv nuclide may be used to stopping power business
if nuclide[:2] == 'Kr':
    dedx_sub_nuclide = 'Kr91'

fraction_escape_up, fraction_escape_down = get_fraction(Shot(134), dedx_sub_nuclide)

print(f"fraction_escape_up: {fraction_escape_up}    fraction_escape_down: {fraction_escape_down}")

fiss_rate_down = np.sum(u238.gamma_induced_fiss_xs.interp(tally_down.energies)*tally_down.dx_per_src*atom_density)
fiss_rate_up = np.sum(u238.gamma_induced_fiss_xs.interp(tally_up.energies)*tally_up.dx_per_src*atom_density)

n_fiss_down = ni_meas_scale * fiss_rate_down
n_fiss_up = ni_meas_scale * fiss_rate_up

if suppress_upstream:
    n_fiss_up *= 0

n_ff_model = (n_fiss_down * fraction_escape_down + n_fiss_up * fraction_escape_up) * ff_yield

print(f'\n\nmodel: {n_ff_model}\nMeasurement: {n_ff_meas}\nefficiency: {100 * n_ff_meas / n_ff_model:.0f}%')
print()
plt.show()