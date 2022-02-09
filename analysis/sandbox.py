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
from JSB_tools import mpl_hist

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

Shot(134).list.plotly(erg_max=2000, remove_baseline=True, convolve_overlay_sigma=1, time_bin_width=40)