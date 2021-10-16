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
eff_220 = 0.1  # efficiency at 218 KeV
#  ======================================================
ni_meas_scale = n_electrons * model_correction
fraction_escape = 0.1
xe_139 = Nuclide.from_symbol("Xe139")

shot = Shot(134)
fit_result = shot.llnl_spe.multi_peak_fit([218.5, 221], fit_window=50)
n_xe_counts = ufloat(fit_result.params['_0amplitude'].value, fit_result.params['_0amplitude'].stderr)
n_xe_meas = n_xe_counts/eff_220/xe_139.decay_gamma_lines[0].intensity

print(fit_result.fit_report())
print(fit_result)
shot.llnl_spe.plot_erg_spectrum()

outp = OutP(Path(__file__).parent/'sims'/'1_inp'/'outp')
tally_up = outp.get_tally('Active up')
tally_down = outp.get_tally('Active down')

atom_density = tally_up.cell.atom_density

u238 = Nuclide.from_symbol('U238')
u238.gamma_induced_fiss_xs.plot()
yields = FissionYields('U238', 'gamma', tally_up.energies)
xe139_yield = np.average(yields.get_yield('Xe139'), weights=tally_up.nominal_fluxes)
print(xe139_yield)

fiss_rate_down = np.sum(u238.gamma_induced_fiss_xs.interp(tally_down.energies)*tally_down.dx_per_src*atom_density)
n_fiss_down = ni_meas_scale * fiss_rate_down * fraction_escape
print('model: ', n_fiss_down)
print('gamma spec: ', n_xe_meas)
plt.show()