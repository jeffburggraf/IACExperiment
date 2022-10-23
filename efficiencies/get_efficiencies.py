import pickle
import matplotlib.pyplot as plt
import numpy as np
from JSB_tools.spe_reader import SPEFile
from pathlib import Path
from JSB_tools.spectra import EfficiencyCalMixin
import uncertainties.unumpy as unp
from JSB_tools.regression import LogPolyFit
from analysis import Shot
from mpant_reader import MPA
from lmfit.model import ModelResult, load_modelresult
from efficiencies.rel_eff import fit_func
# 3 mm vs coffee cup efficiency is a factor of ~1.6
# 0 mm vs 3 mm  efficiency is a factor of ~1.12

pwd = Path(__file__).parent

# with open(Path(__file__).parent/'rel_eff.pickle', 'rb') as f:
rel_model: ModelResult = load_modelresult(str(Path(__file__).parent/'rel_eff.pickle'), {'fit_func': fit_func})

test_shot = Shot(134)
test_spe = test_shot.llnl_spe

with open(Path(__file__).parent/'xy_points.pickle', 'rb') as f:
    xy_points = pickle.load(f)

meas_x = xy_points[0.317]['x']
meas_y = xy_points[0.317]['y']
meas_y_err = unp.std_devs(xy_points[0.317]['y'])
meas_y = unp.nominal_values(meas_y)


angle_x = []
angle_y = []

with open(Path(__file__).parent/'angle_3mm') as f:
    f.readline()
    for line in f.readlines():
        x, _, y = map(float, line.split())
        angle_y.append(y)
        angle_x.append(x)


angle_x = np.array(angle_x)
angle_y = np.array(angle_y)


far_meas_y = unp.nominal_values(xy_points[15.7]['y'])
far_meas_y = far_meas_y*max(meas_y)/max(far_meas_y)

if __name__ == '__main__':
    fig, (ax1, ax2) = plt.subplots(1, 2)

    ax1.plot(angle_x, angle_y, label="Angle Calculation from 15 mm")
    ax1.errorbar(meas_x, meas_y, meas_y_err, label="Measured at 3 mm", ls='None', marker='o')
    ax1.errorbar(xy_points[15.7]['x'],  far_meas_y, label="Measured at 15 mm (normalized)", ls='None', marker='o')
    ax1.legend()

angle_model = LogPolyFit(angle_x, angle_y, 0.05*angle_y, fix_coeffs=[0, 1], order=4)

list_spec = Shot(134).list

base_eff_x = list_spec.erg_centers[np.where(list_spec.erg_centers < 1800)]

base_llnl_effy = unp.nominal_values(angle_model.eval_fit(base_eff_x))

base_iac_effy = base_llnl_effy*unp.uarray(rel_model.eval(x=base_eff_x),
                                          rel_model.eval_uncertainty(x=base_eff_x))


eff_main = EfficiencyCalMixin()
eff_cold = EfficiencyCalMixin()

eff_main.set_efficiency(base_llnl_effy, base_eff_x)
eff_cold.set_efficiency(base_llnl_effy/1.6, base_eff_x)
eff_main.pickle_efficiency(pwd/'eff_main')
eff_cold.pickle_efficiency(pwd/'eff_cold')


nickel_spe_llnl = SPEFile(pwd.parent/'exp_data'/'Nickel'/'Nickel.Spe')
nickel_spe_llnl.set_efficiency(base_llnl_effy*1.1, base_eff_x)
nickel_spe_llnl.pickle_efficiency()


nickel_spe_iac = MPA(pwd.parent/'exp_data'/'Nickel'/'Nickel_iac.mpa')
nickel_spe_iac.set_efficiency(base_iac_effy*1.05, base_eff_x)
nickel_spe_iac.pickle_efficiency()


def set_eff(list_file ,shotnum):
    eff_obj = eff_main

    if 124 <= shotnum <= 134 or 93 <= shotnum <= 101:  # account for 1 cm added distance
        eff_obj = eff_cold

    list_file.path.with_suffix('.eff').unlink(missing_ok=True)

    list_file.eff_path = eff_obj.eff_path



if __name__ == '__main__':
    plt.show()
