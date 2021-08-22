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
from collections import namedtuple

#  =============================================================
path = "their_det/2021-08-17"
reload = True
debug_print = False
false_code = False
pickle_cal = True
#  =============================================================


def cal(rel_data_dir_path, debug=False):
    with open(top_level_data_path/rel_data_dir_path/'peak_points.pickle', 'rb') as f:
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
    y[2] /= (1-(fit.eval_fit(x=121.78)*0.46 + fit.eval_fit(x=867.38)*0.506)).n
    fit = LogPolyFit(x, y, yerr=y_err, order=3, fix_coeffs=[0, 3])

    fit.plot_fit()
    if pickle_cal:
        with open(top_level_data_path/rel_data_dir_path/"cal.pickle", 'wb') as f:
            pickle.dump(fit.fit_result, f)
    return fit.fit_result


if __name__ == '__main__':
    def get_gammas(n: Nuclide, debug=False):
        """
        Return Nuclide.Gammaline instances for each nuclide. Edit this when adding more nuclides.
        Args:
            n:

        Returns:

        """
        n_parent = n
        # add special cases of daughter decays here.
        if n.name == 'Cd109':  # gamma comes from daughter nucleus
            n = Nuclide.from_symbol("Ag109_m1")
        elif n.name == 'Cs137':  # gamma comes from daughter nucleus
            n = Nuclide.from_symbol("Ba137_m1")
        else:
            n = n

        # add special cases of multiple nuclides here
        if n.name == 'Y88':
            cal_gammas = n.decay_gamma_lines[:2]
        elif n.name == 'Eu152':
            cal_gammas = [n.decay_gamma_lines[2], n.decay_gamma_lines[3], n.decay_gamma_lines[7]]
        elif n.name == 'Co57':
            cal_gammas = n.decay_gamma_lines[:1]
        else:
            cal_gammas = n.decay_gamma_lines[:1]

        if debug:
            print(f"Decay gamma lines for {n_parent}:")
            for g in n.decay_gamma_lines:
                print(f"\t{g}")
            print()

        out = namedtuple("gammas", "cal_gammas all_gammas")
        out.cal_gammas = cal_gammas
        out.all_gammas = n.decay_gamma_lines[:]
        return out

    from JSB_tools import Nuclide
    for g in Nuclide.from_symbol('Eu152').decay_gamma_lines:
        print(g)

    def exclude_src_function(s):
        # if 'Eu152' in s:
        #     return False
        return True
    #
    if reload:
        save_eff_points(path, gamma_func=get_gammas, exclude_source_func=exclude_src_function, debug=debug_print)  # Save eff points to file. Change debug to True if needed.
    #

    cal(path, debug_print)
    plt.show()

