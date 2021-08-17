from JSB_tools import Nuclide
import numpy as np
from matplotlib import pyplot as plt
import re
from pathlib import Path
from JSB_tools.list_reader import MaestroListFile
from JSB_tools.spe_reader import SPEFile

data_path = Path.cwd()/'cal_data'/'our_det'/'08_17'

l = MaestroListFile(data_path/'test_beam.Lis')


l.plot_count_rate()
l.plot_sample_ready()
plt.show()
