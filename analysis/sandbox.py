# import re
# import warnings
# import numpy as np
# from JSB_tools.spe_reader import  SPEFile
# from mpant_reader import MPA
# import matplotlib.pyplot as plt
# from pathlib import Path
# from cal_sources import CalSource
# from JSB_tools import Nuclide
# from JSB_tools.nuke_data_tools import GammaLine
# from JSB_tools import mpl_hist
# from JSB_tools.nuke_data_tools.gamma_coince import Levels
# import urllib.request
# from uncertainties import ufloat
# import pandas as pd
# from typing import Tuple
# from uncertainties import unumpy as unp
# from JSB_tools.regression import LogPolyFit
#
# #  CD-109 at 88 KeV
# # spe = SPEFile('/Users/burggraf1/PycharmProjects/IACExperiment/cal_data/our_det/2021-08-17/Cd-109-0.Spe')
# # src= CalSource.get_source('129757')
# # counts = spe.get_counts(87, 90, remove_baseline=True, debug_plot=True)
# # eff = counts / (0.036*src.get_n_decays(spe.livetime, spe.system_start_time))
# # print(sum(eff))
# # spe.plot_erg_spectrum(remove_baseline=True)
# # plt.show()
#
#
# sources = {'Cs137': CalSource.get_source(129792),
#            'Eu152': CalSource.get_source(129753),
#            'Co57': CalSource.get_source('K4-895')
#            }
#
# Ba137_m1 = Nuclide.from_symbol('Ba137_m1')
# co57 = Nuclide.from_symbol('Co57')
# eu152 = Nuclide.from_symbol('Eu152')
#
# #  ===============================================================
# default_window = 3, 5  # first is llnl seconds is iac
# # {'nuclide_name': [(gamma1, window_1), (gamma2, window_2), ...]}
# #     where window_i is either None for default window or a tuple of TWO numbers representing gamma counting window.
# gammas = {'Co57': [(co57.decay_gamma_lines[0], None),
#                    (co57.decay_gamma_lines[1], None)
#                    ],
#           'Cs137': [
#               (Ba137_m1.decay_gamma_lines[0], None)
#                     ],
#           'Eu152': [
#               (eu152.get_gamma_nearest(1085.841), (1084, 1088)),
#               (eu152.get_gamma_nearest(344.279), None),
#               (eu152.get_gamma_nearest(244.698), None)
#                     ]
#           }
# #  ===============================================================
#
# p_data = Path.cwd().parent/'exp_data'
# iac_ni = MPA(p_data/'Nickel'/'Nickel.mpa')
# llnl_ni: SPEFile = SPEFile(p_data/'Nickel'/'Nickel.Spe')
#
# llnl_co57_fri = SPEFile(p_data/'friday'/'EffCalFriday'/'Co57-0.Spe')
# llnl_cs137 = SPEFile(p_data/'friday'/'EffCalFriday'/'Cs137-0.Spe')
# llnl_co57_end_day = SPEFile(p_data/'friday'/'EffCalFriday'/'EndOfDay'/'Co57-0.Spe')
# llnl_eu152 = SPEFile(p_data/'friday'/'EffCalFriday'/'Eu152-0.Spe')
# iac_co57: SPEFile = MPA(p_data/'friday'/'MCA'/'EffCalFriday'/'EndOfDay'/'Co57-0.mpa')
# iac_eu152: SPEFile = MPA(p_data/'friday'/'MCA'/'EffCalFriday'/'Eu152-0.mpa')
# iac_cs137: SPEFile = MPA(p_data/'friday'/'MCA'/'EffCalFriday'/'Cs137-0.mpa')
#
#
# def eff(spe_file: SPEFile):
#     if spe.path.suffix == '.Spe':
#         llnl = True
#     else:
#         llnl = False
#
#     m = re.match("([A-Z][a-z]+[0-9]{1,3}).+", spe_file.path.name)
#     if m:
#         nuclide_name = m.groups()[0]
#     else:
#         assert False
#     src: CalSource = sources[nuclide_name]
#     n_decays = src.get_n_decays(spe_file.livetime, spe_file.system_start_time)
#     x = list(map(lambda x: x[0].erg.n, gammas[nuclide_name]))
#     y = []
#     _default_window = default_window[0] if llnl else default_window[-1]
#
#     for g, window in gammas[nuclide_name]:
#         if window is None:
#             erg_min = g.erg.n-_default_window/2
#             erg_max = g.erg.n+_default_window/2
#         else:
#             assert len(window) == 2
#             erg_min, erg_max = window
#         g: GammaLine
#         counts = spe_file.get_counts(erg_min,erg_max, remove_baseline=True, debug_plot=True)
#         counts = sum(counts)
#         y.append(counts/(n_decays*g.intensity))
#     return nuclide_name, x, y
#
#
# xs, ys = [], []
# plt.figure()
# ax = plt.gca()
#
#
# for spe in [llnl_co57_fri, llnl_cs137, llnl_eu152]:
#     spe.plot_erg_spectrum(remove_baseline=True)
#     name, x, y = eff(spe)
#     ax.errorbar(x, unp.nominal_values(y), unp.std_devs(y), ls='None', marker='o', label=name)
#     xs.extend(x)
#     ys.extend(y)
#
#
# l = Levels('Eu152')
# # ys[4] *= 1.25
# ys = [ufloat(1E-7, 0.15)] + ys
# xs = [35] + xs
# print(ys)
# fit = LogPolyFit(xs, ys, fix_coeffs=[2, 3, 4], order=4)
# fit.plot_fit()
# print(fit.fit_result.fit_report())
# print(fit.eval_fit(778))
# print(l.print_coinc(probability_cut_off=0.01))
#
# plt.show()
#

from scipy.stats.mstats import winsorize
import numpy as np
a = []

print(np.mean(winsorize(a, limits=[0.2,0.2])))