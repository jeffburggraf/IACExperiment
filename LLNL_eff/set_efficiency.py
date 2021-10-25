from JSB_tools.list_reader import MaestroListFile
from JSB_tools.spe_reader import SPEFile, EfficiencyContainer
import pickle
import numpy as np
from analysis import _get_maesto_list_shot_paths
from pathlib import Path


pickle_path = Path(__file__).parent/'3_17mm.pickle'
with open(pickle_path, 'rb') as f:
    llnl_model = pickle.load(f)

llnl_eff = EfficiencyContainer(model=llnl_model)

for shot, p in _get_maesto_list_shot_paths().items():
    list_file = MaestroListFile(p)
    llnl_eff.recalc_effs(list_file.erg_centers)
    list_file.efficiency_meta = llnl_eff
    list_file.pickle()
    print(f"Shot {shot} done")
    # print(p)