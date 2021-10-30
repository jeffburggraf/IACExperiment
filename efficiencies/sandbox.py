import matplotlib.pyplot as plt

from analysis import Shot


s = Shot(120, load_erg_cal=False)
s.list.plotly(remove_baseline=True)
# ax = s.llnl_spe.plot_erg_spectrum(eff_corr=True, make_density=True)
# s.iac_spe.plot_erg_spectrum(eff_corr=True, ax=ax, make_density=True)
# s.iac_spe.plot_efficiency()
#
# plt.show()