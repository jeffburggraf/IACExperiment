import matplotlib.pyplot as plt
from JSB_tools import decay_nuclide
from JSB_tools import Nuclide, FissionYields
from pathlib import Path
from JSB_tools.MCNP_helper.outp_reader import OutP
import scipy
import numpy as np
from openmc.data import Decay

d_jeff = Decay.from_endf('/Users/burggraf1/PycharmProjects/JSB_tools/JSB_tools/nuke_data_tools/endf_files/decay/decay_JEFF/I137')
d_endf = Decay.from_endf('/Users/burggraf1/PycharmProjects/JSB_tools/JSB_tools/nuke_data_tools/endf_files/decay/decay_ENDF/dec-053_I_137.endf')
print()

def f(x):
    ergs = []
    intensities = []
    n = x.spectra['gamma']['discrete_normalization']
    for d in x.spectra['gamma']['discrete']:
        ergs.append(d['energy']*1E-3)
        intensities.append(d['intensity']*n)

    out = list(sorted(zip(ergs, intensities), key=lambda e_i: -e_i[1]))
    return out

print(f(d_jeff))
print(f(d_endf))