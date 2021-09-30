from JSB_tools.MCNP_helper.outp_reader import OutP
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from JSB_tools.MCNP_helper import units
from erg_dist import neutron_ergs, neutron_fluxes
from JSB_tools import mpl_hist

det_area = 6*30*units.inch**2  # in cm^2
foil_area = 1  # cm^2

outp = OutP(Path.cwd()/'0_inp'/'outp')

n_src_neutrons = sum(neutron_fluxes)*foil_area

tally_n = outp.get_cell_by_name('Detector').tallys[0]
tally_p = outp.get_cell_by_name('Detector').tallys[1]
ax = mpl_hist(tally_n.energy_bins, det_area*n_src_neutrons*tally_n.fluxes, label='Neutrons')
mpl_hist(tally_p.energy_bins, det_area*n_src_neutrons*tally_p.fluxes, ax=ax, label="Photons")
ax.set_ylabel("Particles entering perfect detector [Hz]")
ax.set_xlabel("Energy [MeV]")
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.show()