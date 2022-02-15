import numpy as np
from matplotlib import pyplot as plt
from JSB_tools import mpl_hist
from pathlib import Path
from JSB_tools.SRIM import _SRIMConfig, existing_outputs, get_srim_output

from JSB_tools.MCNP_helper import StoppingPowerData
from JSB_tools.MCNP_helper.materials import _IdealGasProperties
#



o = get_srim_output(['He', 'ar'], [1,1], 1E-3, 'xe139', True)

# a =_SRIMConfig(['He', 'ar'], [1,1], 1E-3, 'xe139', True)


print()
for b in _SRIMConfig.all_configs():
    if a == b:
        print(b)
        # print(b.file_name)
# from JSB_tools.SRIM import SRIMTable
#
# gas_density = _IdealGasProperties(['He', 'Ar']).get_density_from_atom_fractions([1, 1], pressure=1.2)
#
#
# srim = SRIMTable(['He', 'Ar'], [1,1], gas_density, 'Xe139', True)
#
# ax = srim.plot_dedx()
#
# mcnp = StoppingPowerData.get_stopping_power("Xe139", material_element_symbols=['Ar', 'He'], density=gas_density,
#                                      material_atom_percents=[1, 1], gas=True, )
#
# plt.show()
