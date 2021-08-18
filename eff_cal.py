"""
Load efficiency data from eff_calibration_setup.py and perform a fit.
Fit can be pickled for later use.

"""
from __future__ import annotations
from matplotlib import pyplot as plt
import pickle
from JSB_tools.regression import LogPolyFit
from eff_calibration_setup import save_eff_points, top_level_data_path
import warnings

#  =============================================================
path = "our_det/08_17"
reload = False
debug_print = True
false_code = True
pickle_cal = True
#  =============================================================


def cal(rel_data_dir_path, debug=False):
    with open(top_level_data_path/rel_data_dir_path/'points.pickle', 'rb') as f:
        _d = pickle.load(f)
        x = _d['x']
        y = _d['y']
        y_err = _d['y_err']
    if debug:
        plt.errorbar(x, y, y_err, ls='None', marker='p')
        plt.xlabel("energy")
        plt.ylabel("Efficiency")
    if false_code:
        y[2] *= 17/20
        warnings.warn("False code in use!")
    fit = LogPolyFit(x, y, yerr=y_err, order=3, fix_coeffs=[0, 3])
    fit.plot_fit()
    if pickle_cal:
        with open(top_level_data_path/rel_data_dir_path/"cal.pickle", 'wb') as f:
            pickle.dump(fit.fit_result, f)
    return fit.fit_result


if __name__ == '__main__':
    if reload:
        save_eff_points(path, False)  # Save eff points to file. Change debug to True if needed.

    cal(path, debug_print)
    plt.show()

