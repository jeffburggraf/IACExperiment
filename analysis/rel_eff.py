"""
Lines to consider using:
974, decays out in time but bad counts
Cs140 at 602 KeV
La-144 @ 397
"""
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


shots_llnl = [85, 86]  # Not including shot 79 since it's short acquisition duration is a bottleneck
shots_iac = [80, 81]
shots_dict: Dict[int, Shot] = {i: Shot(i) for i in shots_iac+shots_llnl}
max_time = min(s.list.total_realtime for s in shots_dict.values())

print("IAC MPA shot times: ", [(num, shots_dict[num].llnl_spe.realtime) for num in shots_iac])
print("IAC Lis shot times: ", [(num, shots_dict[num].list.total_realtime) for num in shots_llnl])

iac_spec = Shot(81).iac_spe
llnl_spec = Shot(85).list.build_spe(0, iac_spec.realtime)

x, y = [], []
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


add_point(218.5, 2, baseline_method='median',  baseline_kwargs={'window_kev': 45})
add_point(397, (1.5, 2.5))
add_point(137.8, 4)
add_point(602.5, 2.5)
add_point(69, 4, baseline_method='median',  baseline_kwargs={'window_kev': 45})
add_point(1428.5, 3)

srt = np.argsort(x)
x = np.array(x)[srt]
y = np.array(y)[srt]


print(x, unp.nominal_values(y))


def fit_func(x, a, b):
    return [a*log(xi*b) for xi in x]


model = Model(fit_func)
params = model.make_params()
params['a'].set(value=0.5)
params['b'].set(value=0.025)

fit = model.fit(x=x, data=unp.nominal_values(y), weights=1.0/unp.std_devs(y), params=params)
# fit.plot_fit()
print(fit.fit_report())


def efficiency_result(x):
    return fit_func(x, ufloat(fit.params['a'].value, fit.params['a'].stderr),
                    ufloat(fit.params['b'].value, fit.params['b'].stderr))


if __name__ == '__main__':
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