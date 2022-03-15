import numpy as np
from matplotlib import pyplot as plt
from JSB_tools import mpl_hist, Nuclide
from pathlib import Path

n = Nuclide.from_symbol('Rb89')
erg = None

if erg is not None:
    es = np.array([x.erg.n for x in n.decay_gamma_lines])
    i = np.argmin(np.abs(erg - es))
    print(f"index: {i}; {n.decay_gamma_lines[i]}\n")

for index, g in enumerate(n.decay_gamma_lines):
    print(index, g)