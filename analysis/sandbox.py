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

data_path = Path(__file__).parent.parent/'exp_data'
bg_iac = MPA(data_path/'tuesday/MCA/BG001.mpa')
bg_llnl = SPEFile(data_path/'tuesday/BG.Spe')

co57_llnl = SPEFile(data_path/'friday'/'EffCalFriday'/'EndOfDay'/'Co57-0.Spe')

ni_iac = MPA(data_path/'Nickel'/'Nickel.mpa')
ni_llnl = SPEFile(data_path/'Nickel'/'Nickel.Spe')

ax_iac = bg_iac.plot_erg_spectrum(make_rate=True, make_density=True)
ax_llnl = bg_llnl.plot_erg_spectrum(make_rate=True, make_density=True)

ni_iac.plot_erg_spectrum(make_rate=True, make_density=True, ax=ax_iac)
ni_llnl.plot_erg_spectrum(make_rate=True, make_density=True, ax=ax_llnl)

co57_llnl.plot_erg_spectrum(make_density=True, make_rate=True, ax=ax_llnl)

plt.show()

