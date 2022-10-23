"""
Lines to consider using:
974, decays out in time but bad counts
Cs140 at 602 KeV
La-144 @ 397
"""
import pickle

import matplotlib.pyplot as plt
import numpy as np
from analysis import Shot, get_merged_time_dependence
# import numpy as np
from JSB_tools import mpl_hist, convolve_gauss
# import matplotlib.pyplot as plt
from JSB_tools.regression import LinearFit
from uncertainties import ufloat
from uncertainties import unumpy as unp
from lmfit.models import LinearModel
from lmfit.models import ExponentialModel
from typing import Dict
from lmfit import Model
from uncertainties.umath import log
from uncertainties import UFloat
from lmfit.model import save_modelresult
#  =======================================
debug = False
#  =======================================


def add_point(erg, window_kev, baseline_method='root', baseline_kwargs=None):
    if isinstance(window_kev, (float, int)):
        window_kev_llnl = window_kev_iac = window_kev
    else:
        assert len(window_kev) == 2
        window_kev_llnl, window_kev_iac = window_kev

    c_iac = iac_spec.get_counts(erg - window_kev_iac, erg + window_kev_iac, remove_baseline=True,
                                baseline_method=baseline_method,
                                debug_plot=debug, baseline_kwargs=baseline_kwargs)

    c_llnl = llnl_spec.get_counts(erg - window_kev_llnl, erg + window_kev_llnl, remove_baseline=True,
                                  baseline_method=baseline_method,
                                  debug_plot=debug, baseline_kwargs=baseline_kwargs)
    x.append(erg)
    y.append(sum(c_iac) / sum(c_llnl))

# def fit_func(x, a, b, c=0):
#     def f(_x):
#         return a*log(_x*b + c * _x**2)
#
#     if isinstance(a, UFloat):
#         return np.array([f(xi) if xi > 0 else 1 for xi in x])
#     else:
#         if x[0] < 0:
#             imin = np.searchsorted(x, 0)
#             out = np.ones_like(x)
#             out[imin:] = f(x[imin:])
#         else:
#             return a*np.log(x*b)
#         return out


def fit_func(x, a, b):
    if isinstance(a, UFloat):
        return np.array([a*log(xi*b) if xi > 0 else 1 for xi in x])  # + xi**2*c
    else:
        if x[0] < 0:
            imin = np.searchsorted(x, 0)
            out = np.ones_like(x)
            out[imin:] = a*np.log(b*_x)
        else:
            return a*np.log(x*b)
        return out


def efficiency_result(x):
    return fit_func(x, ufloat(fit.params['a'].value, fit.params['a'].stderr if fit.params['a'].stderr is not None else 0),
                    ufloat(fit.params['b'].value, fit.params['b'].stderr if fit.params['b'].stderr is not None else 0),
                    )


if __name__ == '__main__':
    # plotly of shots 85 and 86
    # l_85 = Shot(85).list
    # l_85 += Shot(86).list
    # l_85.plotly(time_bin_width=20)
    shots_llnl = [85, 86]  # Not including shot 79 since it's short acquisition duration is a bottleneck
    shots_iac = [80, 81]
    shots_dict: Dict[int, Shot] = {i: Shot(i) for i in shots_iac + shots_llnl}
    max_time = min(s.list.total_realtime for s in shots_dict.values())

    print("Using max time: ", max_time)
    print("IAC MPA shot times: ", [(num, shots_dict[num].llnl_spe.realtime) for num in shots_iac])
    print("IAC Lis shot times: ", [(num, shots_dict[num].list.total_realtime) for num in shots_llnl])

    ax = Shot(80).iac_spe.plot_erg_spectrum(make_rate=True, leg_label='Shot 80', eff_corr=False)
    Shot(81).iac_spe.plot_erg_spectrum(ax=ax, make_rate=True, leg_label='Shot 81', eff_corr=False)
    ax.set_title("Two equiv. shots. (IAC det.)")

    iac_spec = Shot(81).iac_spe
    llnl_spec = Shot(85).list.build_spe(0, iac_spec.realtime)

    x, y = [], []

    add_point(218.5, 2, baseline_method='median', baseline_kwargs={'window_kev': 45})
    add_point(397, (1.5, 2.5))
    add_point(137.8, 4)
    add_point(602.5, 2.5)
    add_point(974.5, 2.5, baseline_method='median', baseline_kwargs={'window_kev': 50})
    add_point(697, (2, 3.3))
    add_point(69, 3, baseline_method='median', baseline_kwargs={'window_kev': 45})
    add_point(1428.5, 3, baseline_method='median', baseline_kwargs={'window_kev': 45})

    srt = np.argsort(x)
    x = np.array(x)[srt]
    y = np.array(y)[srt]

    model = Model(fit_func)
    params = model.make_params()
    params['a'].set(value=0.5)
    params['b'].set(value=0.025)

    fit = model.fit(x=x, data=unp.nominal_values(y), weights=1.0 / unp.std_devs(y), params=params)
    # fit.plot_fit()

    # with open("rel_eff.pickle", 'wb') as f:
    #     pickle.dump(fit, f)
    save_modelresult(fit, 'rel_eff.pickle')

    print(fit.fit_report())

    plt.figure()
    plt.errorbar(x, unp.nominal_values(y), unp.std_devs(y), marker='o', ls='None')
    plt.xlabel("Energy [KeV]")
    plt.ylabel("IAC/LLNL")
    _x = np.arange(40, 2000, 10)
    plt.plot(_x, unp.nominal_values(efficiency_result(_x)), label='Fit')
    plt.legend()

    ax = iac_spec.plot_erg_spectrum(make_density=True, leg_label='IAC', remove_baseline=True)
    llnl_spec.plot_erg_spectrum(make_density=True, leg_label='LLNL', ax=ax, remove_baseline=True)
    ax.legend()
    # iac_spec.set_efficiency(np.interp(iac_spec.energies, x, unp.nominal_values(y)))
    plt.show()