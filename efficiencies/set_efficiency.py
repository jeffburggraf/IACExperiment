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
from rel_eff import fit_func
from lmfit.model import ModelResult


# / ===========
pickle_llnl = True
pickle_iac = True
#  ==============================
pickle_path = Path(__file__).parent/'3_17mm.pickle'
with open(pickle_path, 'rb') as f:
    llnl_model: ModelResult = pickle.load(f)

with open( Path(__file__).parent/'rel_eff.pickle', 'rb') as f:
    rel_model = pickle.load(f)


def get_iac_effs(ergs):
    return llnl_model.eval(x=ergs) * unp.uarray(rel_model.eval(x=ergs), rel_model.eval_uncertainty(x=ergs))


iac_shots = _get_mpant_mca_shot_paths()

nickel_path = Path(__file__).parent.parent/'exp_data'/'Nickel'
llnl_nickel = SPEFile(nickel_path/'Nickel.Spe')
llnl_nickel.set_useful_energy_range(40)
llnl_nickel.eff_model = llnl_model
llnl_nickel.pickle_eff()

iac_nickel = MPA(nickel_path/'Nickel.mpa')
iac_nickel.set_useful_energy_range(40)
iac_nickel.effs = get_iac_effs(iac_nickel.energies)
iac_nickel.pickle_eff()

for shot, p in _get_maesto_list_shot_paths().items():
    # if shot < 119:
    #     continue
    print(shot)

    # if pickle_llnl:
    #     list_file = MaestroListFile(p)
    #     list_file.set_useful_energy_range(40)
    #     list_file.eff_model = llnl_model
    #     if 124 <= shot <= 134:  # account for 1 cm added distance
    #         list_file.eff_scale = 0.86
    # else:
    list_file = MaestroListFile.from_pickle(p, load_erg_cal=False)
    list_file.set_useful_energy_range(40)
    list_file.eff_model = llnl_model
    if 124 <= shot <= 134:  # account for 1 cm added distance
        list_file.eff_scale = 0.86

    # if pickle_llnl:
    #     list_file.pickle()
    list_file.pickle_eff()

    if shot in iac_shots:
        iac_spe = MPA(iac_shots[shot])
        iac_spe.set_useful_energy_range(40)
        iac_spe.effs = list_file.interp_eff(iac_spe.energies)*\
                       unp.uarray(rel_model.eval(x=iac_spe.energies), rel_model.eval_uncertainty(x=iac_spe.energies))
        iac_spe.pickle_eff()
        iac_spe.pickle()
        # plt.show()
    else:
        continue
    #  MaestroListFile.from_pickle(p)
