import numpy as np
from JSB_tools.spe_reader import SPEFile
from matplotlib import pyplot as plt


s = SPEFile('/Users/burggraf1/PycharmProjects/IACExperiment/cal_data/our_det/08_17/Eu-152-0.Spe')

s.get_spectrum_hist().plot()
plt.show()