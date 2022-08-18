import matplotlib.pyplot as plt
from JSB_tools.maestro_reader import MaestroListFile
from JSB_tools.spe_reader import SPEFile
import timeit
import cProfile
from lmfit.models import GaussianModel, LinearModel
from lmfit.model import CompositeModel
import numpy as np
from uncertainties import ufloat, nominal_value
import pickle
from JSB_tools.nuke_data_tools import FissionYields, DecayNuclide
from JSB_tools import Nuclide, mpl_hist_from_data, flatten
#!/usr/bin/env python
from JSB_tools.MCNP_helper.outp_reader import OutP
from pathlib import Path

x = np.linspace(0, 10, 100)
y = x**2

# fig, ax = plt.subplots(1, 1)
# ax.plot(x, y)
#
# with open("delete", 'wb') as f:
#     pickle.dump((fig, ax), f)

with open("delete", 'rb') as f:
    (fig, ax) = pickle.load(f)

plt.show()




