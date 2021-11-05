import matplotlib.pyplot as plt

from analysis import Shot


s = Shot(120, load_erg_cal=True)
s.iac_spe.unpickle_eff()
ax = s.iac_spe.plot_efficiency()
s.llnl_spe.plot_efficiency(ax=ax, label="LLNL")
plt.show()