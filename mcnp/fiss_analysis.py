import warnings

from JSB_tools.nuke_data_tools import Nuclide
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

c_per_second = (192/3.0)*1E-6
charge_per_electron = 1.602E-19
n_electrons = 3*c_per_second / charge_per_electron

#  ======================================================
model_correction = 1.0/3.15
eff_220 = 0.07  # efficiency at 218 KeV. Difference due to increased distance is expected to be 75% via solid angle.
#  ======================================================
ni_meas_scale = n_electrons * model_correction
fraction_escape = 0.1
xe_139 = Nuclide.from_symbol("Xe139")

shot = Shot(131)
fit_result = shot.llnl_spe.multi_peak_fit([219, 221], fit_window=30, debug_plot=True)
print(fit_result.fit_report())

amp_error = fit_result.params['_0amplitude'].stderr
if amp_error is None:
    amp_error = fit_result.params['_0amplitude'].value*0.1
    warnings.warn("Fit failed to calculate errors!")

n_xe_counts = ufloat(fit_result.params['_0amplitude'].value, amp_error)
n_xe_meas = n_xe_counts/eff_220/xe_139.decay_gamma_lines[0].intensity


shot.llnl_spe.plot_erg_spectrum()

outp = OutP(Path(__file__).parent/'sims'/'1_inp'/'outp')
tally_up = outp.get_tally('Active up')
tally_down = outp.get_tally('Active down')
tally_down.plot()
tally_up.plot()


atom_density = tally_up.cell.atom_density

u238 = Nuclide.from_symbol('U238')
u238.gamma_induced_fiss_xs.plot()
yields = FissionYields('U238', 'gamma', tally_down.energies)
yields.plot('Xe139',)
xe139_yield = np.average(yields.get_yield('Xe139'), weights=tally_down.nominal_fluxes)

tot = None
for y in yields.yields.values():
    if tot is None:
        tot = y.copy()
    else:
        tot += y
print(tot)

plt.plot(tally_down.energies, u238.gamma_induced_fiss_xs.interp(tally_down.energies))

fiss_rate_down = np.sum(u238.gamma_induced_fiss_xs.interp(tally_down.energies)*tally_down.dx_per_src*atom_density)
n_fiss_down = ni_meas_scale * fiss_rate_down * fraction_escape
print('model: ', n_fiss_down)
print('gamma spec: ', n_xe_meas)
plt.show()