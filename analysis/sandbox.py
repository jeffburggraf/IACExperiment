# import re
# import warnings
# import numpy as np
# from JSB_tools.spe_reader import  SPEFile
from mpant_reader import MPA
import matplotlib.pyplot as plt
from pathlib import Path
import datetime
import numpy as np
from JSB_tools.spe_reader import SPEFile
from cal_sources import CalSource
from JSB_tools.nuke_data_tools.gamma_coince import Levels
from JSB_tools.nuke_data_tools import FissionYields
from analysis import Shot
from JSB_tools.nuke_data_tools.gamma_spec import gamma_search
#
# yields = FissionYields('U238', 'gamma', [5], )
#
#
# def weight(n):
#     return np.mean(yields.get_yield(n))
#
#
# for g in gamma_search(296.5, 2, 15, 400, nuclide_weighting_function=weight):
#     print(g)


ax = None
for i in range(129, 135):
    shot = Shot(i)
    label = shot.__repr__(['he_flow', 'ar_flow', 'shotnum'])
    kwargs = {'time_max': 310, 'time_min': 0, 'remove_baseline': True, 'erg_min': 200, 'erg_max': 250}
    if ax is None:
        ax =shot.list.plot_erg_spectrum(**kwargs, label=label)
        continue

    shot.list.plot_erg_spectrum(ax=ax, **kwargs, label=label)

# Shot(i).list.plotly(remove_background=True)

plt.show()