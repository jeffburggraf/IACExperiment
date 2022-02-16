import numpy as np
from matplotlib import pyplot as plt
from JSB_tools import mpl_hist
from pathlib import Path
from JSB_tools.SRIM import _SRIMConfig, existing_outputs, find_SRIM_run
from JSB_tools.MCNP_helper import StoppingPowerData
from JSB_tools.MCNP_helper.materials import _IdealGasProperties
from JSB_tools.SRIM import find_SRIM_run


gas_density = _IdealGasProperties(['He', 'Ar']).get_density_from_atom_fractions([1, 1], pressure=1.2)


srim = find_SRIM_run(['He', 'Ar'], [1,1], gas_density, 'Xe139', True)

ax = srim.plot_dedx()

mcnp = StoppingPowerData.get_stopping_power("Xe139", material_element_symbols=['Ar', 'He'], density=gas_density,
                                     material_atom_percents=[1, 1], gas=True, )

mcnp.plot_dedx(ax=ax, label='MCNP')
ax.legend()
plt.show()
