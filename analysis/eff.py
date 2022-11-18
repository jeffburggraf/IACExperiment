import numpy as np

from JSB_tools.spe_reader import SPEFile
from mpant_reader import MPA
from matplotlib import pyplot as plt
from cal_sources import CalSource
from datetime import datetime
from JSB_tools.nuke_data_tools import Nuclide
from uncertainties import ufloat
from lmfit.models import ExponentialModel
from uncertainties import unumpy as unp
from JSB_tools import mpl_hist
from JSB_tools.maestro_reader import MaestroListFile
from pathlib import Path


def sanbox1(llnl=True):
    if not llnl:
        base_dir = Path("MCA") / "EffCalFriday"
    else:
        base_dir = "EffCalFriday"

    print(Nuclide('Ni57').decay_gamma_lines)
    co57_src = CalSource.get_source('K4-896')
    cs137_src = CalSource.get_source('129792')
    y88_src = CalSource.get_source('190607000')

    co_57_g1 = Nuclide('Co57').decay_gamma_lines[0]
    co_57_g2 = Nuclide('Co57').decay_gamma_lines[1]
    y88_g1 = Nuclide('Y88').decay_gamma_lines[0]
    ba_137m_g1 = Nuclide('Ba137_m1').decay_gamma_lines[0]
    ni_g1377, ni_g127 = Nuclide('Ni57').decay_gamma_lines[:2]

    print(ba_137m_g1)

    p = Path('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/')

    cls = MPA if not llnl else SPEFile
    s = '.mpa' if not llnl else '.Spe'
    # llnl_co57: SPEFile = cls((p/base_dir/'Co57-0').with_suffix(s))
    llnl_co57: SPEFile = cls((p/base_dir/'EndOfDay/Co57-0').with_suffix(s))
    llnl_cs137: SPEFile = cls((p/base_dir/'Cs137-0').with_suffix(s))
    llnl_eu152: SPEFile = cls((p/base_dir/'Eu152-0').with_suffix(s))
    llnl_y88: SPEFile = cls((p/base_dir/'Y88-0').with_suffix(s))
    if llnl:
        llnl_ni: SPEFile = cls((p.parent/'Nickel'/'Nickel').with_suffix(s))
    else:
        llnl_ni: SPEFile = cls((p.parent/'Nickel'/'Nickel').with_suffix(s))

    llnl_ni.plot_erg_spectrum(make_rate=True)
    llnl_cs137.plot_erg_spectrum(make_rate=True)

    counts1 = sum(llnl_co57.get_counts(co_57_g1.erg.n - 1.5, co_57_g1.erg.n + 1.5, remove_baseline=True,
                                       debug_plot=True, deadtime_corr=False))
    counts2 = sum(llnl_co57.get_counts(co_57_g2.erg.n - 1.5, co_57_g2.erg.n + 1.5, remove_baseline=True,
                                       debug_plot=True, deadtime_corr=False))
    counts3 = sum(llnl_cs137.get_counts(ba_137m_g1.erg.n - 3.5, ba_137m_g1.erg.n + 2.5, remove_baseline=True,
                                        debug_plot=True, deadtime_corr=False))
    counts4 = sum(llnl_y88.get_counts(y88_g1.erg.n - 4.5, y88_g1.erg.n + 4.5, remove_baseline=True, debug_plot=True,
                                      deadtime_corr=False))

    counts_ni_127 = sum(llnl_ni.get_counts(ni_g127.erg.n - 3.5, ni_g127.erg.n + 3.5, remove_baseline=True,
                                           debug_plot=True, baseline_method='median',
                                           baseline_kwargs={'window_kev': 50}, deadtime_corr=True))
    counts_ni_1377 = sum(llnl_ni.get_counts(ni_g1377.erg.n - 3.5, ni_g1377.erg.n + 3.5, remove_baseline=True,
                                            debug_plot=True, baseline_method='median', deadtime_corr=False))

    eff1 = counts1/(co57_src.get_n_decays(llnl_co57.livetime) * co_57_g1.intensity)
    eff2 = counts2/(co57_src.get_n_decays(llnl_co57.livetime) * co_57_g2.intensity)
    eff3 = counts3/(cs137_src.get_n_decays(llnl_cs137.livetime) * ba_137m_g1.intensity)
    # eff4 = counts4/(y88_src.get_n_decays(llnl_y88.livetime) * y88_g1.intensity)
    ni_rel_eff = (counts_ni_1377 / ni_g1377.intensity) / (counts_ni_127 / ni_g127.intensity)

    eff_ni_127 = ufloat(np.interp(127.16, [co_57_g1.erg.n, co_57_g2.erg.n], [eff1.n, eff2.n]),
                        np.interp(127.16, [co_57_g1.erg.n, co_57_g2.erg.n], [eff1.std_dev, eff2.std_dev]))
    eff_ni_1300 = ni_rel_eff * eff_ni_127

    print(f'Ni eff at 127.16: {eff_ni_127}')

    x = [co_57_g1.erg.n, co_57_g2.erg.n, ba_137m_g1.erg.n, ni_g1377.erg.n]
    y = [eff1, eff2, eff3, eff_ni_1300]

    yerr = unp.std_devs(y)
    y = unp.nominal_values(y)

    model = ExponentialModel()
    params = model.guess(x=x, data=y)
    fit = model.fit(data=y, x=x, weights=1.0/yerr, params=params)

    plt.figure()
    _x = np.linspace(100, 2000, 10000)
    plt.plot(_x, fit.eval(x=_x, params=params))
    plt.plot(x, y, ls='None', marker='o')
    print(fit.fit_report())
    eu152 = Nuclide('Eu152')
    ax_eu = eu152.plot_decay_gamma_spectrum(min_intensity=0.02)
    ax_y88 = Nuclide('Y88').plot_decay_gamma_spectrum(min_intensity=0.02)

    eu_x = []
    eu_y = []

    y88_x = []
    y88_y = []
    for g in eu152.decay_gamma_lines:
        if g.intensity > 0.02:
            c = llnl_eu152.get_counts(g.erg.n-3, g.erg.n+3, remove_baseline=False, deadtime_corr=True, debug_plot=False)
            c -= llnl_eu152.get_baseline_median(g.erg.n - 3, g.erg.n + 3, window_kev=90)
            # plt.figure()
            # plt.plot(unp.nominal_values(c))
            # plt.title(f"{g.erg}")
            c = sum(c)
            print(g.erg, fit.eval(x=g.erg.n, params=params), params)
            c = c/fit.eval(x=g.erg.n, params=params)
            eu_x.append(g.erg.n)
            eu_y.append(c)

    for g in Nuclide('Y88').decay_gamma_lines:
        if g.intensity > 0.02:
            c = llnl_y88.get_counts(g.erg.n-3, g.erg.n+3, remove_baseline=False, deadtime_corr=True, debug_plot=False)
            c -= llnl_y88.get_baseline_median(g.erg.n - 3, g.erg.n + 3, window_kev=90)
            # plt.figure()
            # plt.plot(unp.nominal_values(c))
            # plt.title(f"{g.erg}")
            c = sum(c)
            print(g.erg, fit.eval(x=g.erg.n, params=params), params)
            c = c/fit.eval(x=g.erg.n, params=params)
            y88_x.append(g.erg.n)
            y88_y.append(c)

    eu_y = unp.nominal_values(eu_y)
    eu_y *= sum(filter(lambda x: x > 0.02, map(lambda x: x.intensity.n, eu152.decay_gamma_lines)))/sum(eu_y)
    ax_eu.plot(eu_x, eu_y, ls='None', marker='o')

    y88_y = unp.nominal_values(y88_x)
    y88_y *= sum(filter(lambda x: x > 0.02, map(lambda x: x.intensity.n, Nuclide('Y88').decay_gamma_lines)))/sum(y88_y)
    ax_y88.plot(y88_x, y88_y, ls='None', marker='d')

    plt.show()


sanbox1(llnl=False)
# mpa_nickel = MPA('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/MCA/Nickel.mpa')
# mpa_bg = MPA('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/tuesday/MCA/BG001.mpa')
#
# counts, bins = mpa_bg.get_counts(return_bin_edges=True, baseline_method='median',
#                                      baseline_kwargs={'window_kev': 10}, remove_baseline=True, make_rate=True,)
# ax = mpl_hist(bins, counts, label='bg')
#
# counts, bins = mpa_nickel.get_counts(return_bin_edges=True, baseline_method='median',
#                                      baseline_kwargs={'window_kev': 10}, remove_baseline=True, make_rate=True)
# mpl_hist(bins, counts, ax=ax, alpha=0.7)
#
# iac_nickel_lis = MaestroListFile('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/Nickel/Nickel.Lis')
# iac_nickel_lis.pickle()
#
# iac_nickel_lis.plotly(time_step=1000, time_bin_width=500)
#
#
#
# plt.show()