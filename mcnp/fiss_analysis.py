"""
Thoughts Nov 1:
    Is Kr89 really Kr89? Check with better statistics by summiong cold filter shots
    How to do decay correction?
    Make more options to disable certain corrections.

"""
import warnings
from JSB_tools.nuke_data_tools import Nuclide, DecayedFissionYields, DecayNuclide
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

gas_mats = {134: {'gasses': ['Ar'], 'fractions': [1], 'pressure': 1.4},
            131: {'gasses': ['He'], 'fractions': [0.5], 'pressure': 1.4}}
#  ======================================================
model_correction = 0.5*(ufloat(1.51, 0.05) + ufloat(1.04, 0.05))  #  from Ni activation model_avg/meas (average of close and far Ni)
# model_correction = 1.0/1.5
nuclide = 'Sr94'  # Kr89, Sb132, Sr94
shot_num = 134
suppress_upstream = True  # Do or don't include FF's which escape from upstream of foil
do_decay_corr = False
assume_trans_time = 20
gamma_index = 0
#  ======================================================
gas_args = gas_mats[shot_num]
shot = Shot(shot_num)
shot.list.plot_sample_ready()

# shot.list.plotly(remove_baseline=True)
shot.list.plot_erg_spectrum(eff_corr=True)

ni_meas_scale = n_electrons / model_correction

ff = Nuclide(nuclide)


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

outp = OutP(Path(__file__).parent/'sims'/'du_shot131'/'outp')
tally_up = outp.get_f4_tally('Active up')
tally_down = outp.get_f4_tally('Active down')

u238 = Nuclide('U238')

neg_weight = sum(u238.gamma_induced_fiss_xs.interp(tally_up.energies)*tally_up.dx_per_src)
pos_weight = sum(u238.gamma_induced_fiss_xs.interp(tally_down.energies)*tally_down.dx_per_src)
print("pos/neg u foil fission rates: ", pos_weight, neg_weight)

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
print("Model correction: ", model_correction)

amp_error = fit_result.params['_0amplitude'].stderr
if amp_error is None:
    amp_error = fit_result.params['_0amplitude'].value*0.1
    warnings.warn("Fit failed to calculate errors!")

n_ff_meas = ufloat(fit_result.params['_0amplitude'].value, amp_error) / gamma_line.intensity
n_ff_meas *= decay_corr

# shot.llnl_spe.plot_erg_spectrum()




atom_density = tally_up.cell.atom_density

# u238.gamma_induced_fiss_xs.plot()
yields = FissionYields('U238', 'gamma', tally_down.energies, independent_bool=True)
# yields.plot(nuclide)


ff_yield = np.average(yields.get_yield(nuclide),
                      weights=tally_down.nominal_fluxes*u238.gamma_induced_fiss_xs.interp(tally_down.energies))


print('Atom density: ', atom_density)
print(f'{nuclide} independent yield: ', ff_yield)


dedx_sub_nuclide = nuclide  # in case an equiv nuclide may be used to stopping power business
if nuclide[:2] == 'Kr':
    dedx_sub_nuclide = 'Kr91'

fraction_escape_up, fraction_escape_down = get_fraction(Shot(shot_num), dedx_sub_nuclide)

print(f"fraction_escape_up: {fraction_escape_up}    fraction_escape_down: {fraction_escape_down}")

u238.gamma_induced_fiss_xs.plot()
fiss_rate_down = np.sum(u238.gamma_induced_fiss_xs.interp(tally_down.energies)*tally_down.dx_per_src*atom_density)
fiss_rate_up = np.sum(u238.gamma_induced_fiss_xs.interp(tally_up.energies)*tally_up.dx_per_src*atom_density)

n_fiss_down = ni_meas_scale * fiss_rate_down
n_fiss_up = ni_meas_scale * fiss_rate_up

time_bins = np.arange(0, 300, 5)
shot.list.plot_time_dependence(energy=fit_ergs[0], bins=time_bins, signal_window_kev=6, make_rate=True, eff_corr=True)

if suppress_upstream:
    n_fiss_up *= 0

n_ff_model = (n_fiss_down * fraction_escape_down + n_fiss_up * fraction_escape_up) * ff_yield

print(f'\n\nmodel: {n_ff_model}\nMeasurement: {n_ff_meas}\nefficiency: {100 * n_ff_meas / n_ff_model:.0f}%')
print()
plt.show()