import numpy as np
from JSB_tools.spe_reader import SPEFile
from matplotlib import pyplot as plt
from JSB_tools import Nuclide

n = Nuclide.from_symbol('Eu152')
print(n.decay_daughters)
for g in n.decay_gamma_lines:
    print(g)
# s = SPEFile('/Users/burggraf1/PycharmProjects/IACExperiment/cal_data/our_det/2021-08-17/Eu-152-0.Spe')
#
# s.get_spectrum_hist().plot()
# plt.show()