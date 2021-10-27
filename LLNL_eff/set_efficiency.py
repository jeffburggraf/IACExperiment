import lmfit.model
from analysis import MPA
from JSB_tools.list_reader import MaestroListFile
from JSB_tools.spe_reader import SPEFile
import pickle
import numpy as np
from analysis import _get_maesto_list_shot_paths, _get_mpant_mca_shot_paths
from pathlib import Path
from uncertainties import unumpy as unp
from matplotlib import pyplot as plt
from analysis.rel_eff import fit as iac_rel_fit


# / ===========
pickle_llnl = False
pickle_iac = True
#  ==============================
pickle_path = Path(__file__).parent/'3_17mm.pickle'
with open(pickle_path, 'rb') as f:
    llnl_model = pickle.load(f)


iac_shots = _get_mpant_mca_shot_paths()

for shot, p in _get_maesto_list_shot_paths().items():
    if pickle_llnl:
        list_file = MaestroListFile(p)
        list_file.set_useful_energy_range(40)
        list_file.eff_model = llnl_model
    else:
        list_file = MaestroListFile.from_pickle(p)

    if pickle_llnl:
        list_file.pickle()

    if shot in iac_shots:
        iac_spe = MPA(iac_shots[shot])
        iac_spe.effs = list_file.interp_eff(iac_spe.energies)*\
                       unp.uarray(iac_rel_fit.eval(x=iac_spe.energies), iac_rel_fit.eval_uncertainty(x=iac_spe.energies))
    else:
        continue
    print()
    #  MaestroListFile.from_pickle(p)
