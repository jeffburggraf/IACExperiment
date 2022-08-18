import matplotlib.pyplot as plt
from JSB_tools import DecayNuclide
from JSB_tools import Nuclide, FissionYields
from pathlib import Path
from JSB_tools.MCNP_helper.outp_reader import OutP
import scipy
import numpy as np
from openmc.data import Decay
from JSB_tools.nuke_data_tools.nudel import LevelScheme

n = LevelScheme("Co57" )

print()