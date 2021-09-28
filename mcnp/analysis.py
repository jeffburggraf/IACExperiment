from JSB_tools.nuke_data_tools import Nuclide
from matplotlib import pyplot as plt
from pathlib import Path
from JSB_tools.MCNP_helper.outp_reader import OutP
import numpy as np

c_per_second = 192/3.0*1E-6
charge_per_electron = 1.602E-19
n_electrons = 3*c_per_second / charge_per_electron


ni58 = Nuclide.from_symbol('Ni58')
ni57 = Nuclide.from_symbol("Ni57")
u238 = Nuclide.from_symbol('U238')

ni_xs = ni58.get_incident_gamma_daughters()['Ni57'].xs
ni_xs.plot()

out = OutP(Path.cwd()/'sims'/'0_inp'/'outp')
tally_chamber = out.get_cell_by_name('Chamber target').get_tally()

vcd_ni_tally = out.get_cell_by_name('VCD nickel').get_tally()

tally_chamber.plot()
vcd_ni_tally.plot()

ni_per_src_chamber = np.sum(tally_chamber.cell.atom_density*ni_xs.interp(tally_chamber.energies)*tally_chamber.total_dx_per_src)
tot_ni = n_electrons*np.sum(tally_chamber.cell.atom_density*ni_xs.interp(tally_chamber.energies)*tally_chamber.total_dx_per_src)
print(tot_ni)
plt.show()