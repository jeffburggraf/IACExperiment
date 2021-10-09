"""
Analysis of multiple filters.

Notes:
    Issues with detector warmup continued. This makes absolute comparison difficult.
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

shots = map(Shot, [38, 137, 138, 139])


x = []
y = []


for shot in shots:
    window = 216.5, 221.5
    shot.llnl_spe.plot_erg_spectrum(remove_baseline=True)
    spe_llnl = shot.list.build_spe(0, shot.iac_spe.realtime)
    counts_llnl, bins_llnl = spe_llnl.get_counts(*window, remove_baseline=True,
                                                  return_bin_edges=True, make_rate=True)
    llnl_b_widths = bins_llnl[1:] - bins_llnl[:-1]

    counts_iac, bins_iac = shot.iac_spe.get_counts(*window, remove_baseline=True,
                                                   return_bin_edges=True, make_rate=True)
    iac_b_widths = bins_iac[1:] - bins_iac[:-1]

    ax = mpl_hist(bins_llnl, counts_llnl/llnl_b_widths, label='LLNL')
    mpl_hist(bins_iac, counts_iac/iac_b_widths, label=f'IAC - {shot.filters_before_iac_det} filters', ax=ax)
    ratio = sum(counts_iac)/sum(counts_llnl)
    print(shot.filters_before_iac_det, ratio)
    if shot.shotnum >=137:
        x.append(shot.filters_before_iac_det)
        y.append(ratio)
    ax.set_title(f'Shot{shot.shotnum}')
    ax_time = None


plt.figure()
plt.errorbar(x, unp.nominal_values(y), unp.std_devs(y), ls='None', marker='o')

model = ExponentialModel()
params = model.guess(data=unp.nominal_values(y), x=x)
fit_result = model.fit(unp.nominal_values(y), params=params, x=x, weights=1.0/unp.std_devs(y))
print(fit_result.fit_report())

_x = np.linspace(1, 2.3*params['decay'].value, 100)

plt.plot(_x, fit_result.eval(params=params, x=_x), label='Fit')

plt.legend()
plt.xticks(np.arange(1, int(max(_x)+1)))
plt.ylabel("(counts downstream)/(counts upstream)")
plt.xlabel("# filters")
plt.show()
