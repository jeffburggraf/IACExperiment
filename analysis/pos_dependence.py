import numpy as np
import matplotlib.pyplot as plt

from JSB_tools.spe_reader import SPEFile

from cal_sources import CalSource
from JSB_tools.nuke_data_tools.gamma_coince import Levels
l = Levels('Ni57')

n = 'Co57'
erg = {'Co57': 122, 'Cs137': 661.8, 'Eu152':121}[n]
dir = "/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/wednesday/EffCalWednesday/"
s_0 = SPEFile(f'{dir}{n}-0.Spe')
s_1 = SPEFile(f'{dir}{n}-1.Spe')
s_2 = SPEFile(f'{dir}{n}-2.Spe')
ax = s_0.plot_erg_spectrum(leg_label='0', make_rate=True)
s_1.plot_erg_spectrum(leg_label='1', ax=ax, make_rate=True)
s_2.plot_erg_spectrum(leg_label='2', ax=ax,  make_rate=True)

c_0 = sum(s_0.get_counts(erg-2, erg+2, debug_plot=True, make_rate=True))
c_1 = sum(s_1.get_counts(erg-2, erg+2, make_rate=True))
c_2 = sum(s_2.get_counts(erg-2, erg+2, make_rate=True))
factor = 1.0/(np.mean([c_2, c_1, c_0])/c_0)
print(factor)
plt.show()