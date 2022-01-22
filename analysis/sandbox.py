# import re
# import warnings
# import numpy as np
# from JSB_tools.spe_reader import  SPEFile
from mpant_reader import MPA
import matplotlib.pyplot as plt
from pathlib import Path
import datetime
import numpy as np
from JSB_tools import mpl_hist
from JSB_tools.spe_reader import SPEFile
from cal_sources import CalSource
from JSB_tools.nuke_data_tools.gamma_coince import Levels
from JSB_tools.nuke_data_tools import FissionYields
from analysis import Shot
from JSB_tools.nuke_data_tools.gamma_spec import gamma_search
from JSB_tools import Nuclide


# l = None
#
#
# for s in Shot.find_shots(flow='100001', he_flow=0.25, ar_flow=0.25, foil_pos='upstream', cold_filter=False):
#     if l is None:
#         l = s.list
#     else:
#         l += s.list
#     print(s)


# l.plotly(erg_min=50)
# l.plot_time_dependence(1413.3, np.arange(0, 500, 5), signal_window_kev=4.5)

l = None
for s in Shot.find_shots(cold_filter=True):
    if l is None:
        l = s.list
    else:
        l += s.list
    print(s)

l.plot_erg_spectrum(erg_min=60, time_max=320, eff_corr=True, remove_baseline=True)
# l.plotly(erg_min=60, erg_max=2300, time_bin_width=7)
# l.plot_time_dependence(1413.3, np.arange(0, 500, 5), signal_window_kev=4.5)

plt.show()
