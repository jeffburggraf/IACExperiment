"""
For 0.5 total gas flow with pattern 010010, the chamber clearing time is about 10 seconds
"""
import matplotlib.pyplot as plt
import numpy as np
from analysis import Shot, get_merged_time_dependence
# import numpy as np
from JSB_tools import mpl_hist, convolve_gauss
# import matplotlib.pyplot as plt
from JSB_tools.regression import LinearFit
from uncertainties import ufloat
import uncertainties.unumpy as unp
from lmfit.models import LinearModel

#  =============================================
shots_dict = {0.25: [9, 45, 26, 27, 22, 23, 122], 0.5: [24, 25, 28, 29, 30, 31]}
shots = shots_dict[0.25]
b_width = 7
flow = 0.25
erg = 805
window = 4
#  =============================================

shots_dict = {k: list(map(Shot, value)) for k, value in shots_dict.items()}
data = {}

max_time = max(map(lambda x: x.max_time, shots_dict[0.25]))
time_bins = np.arange(0, max_time, b_width)

for shot in shots_dict[flow]:
    shot: Shot
    time_dep, _, _ = shot.list.get_time_dependence(erg, debug_plot=False, bins=time_bins, signal_window_kev=window,
                                                   bg_window_kev=13)
    try:
        n = data[shot.tube_len]["n"]
        data[shot.tube_len]["values"] = (data[shot.tube_len]["values"]*n + time_dep)/(n + 1)
        data[shot.tube_len]["n"] += 1
    except KeyError:
        data[shot.tube_len] = {'values': time_dep, 'n': 1}

for k in data.keys():
    s = sum(data[k]['values'])
    data[k]['values'] /= s

plt.figure()

for length, data in data.items():
    mpl_hist(time_bins, unp.nominal_values(data['values']), ax=plt.gca(), label=length, poisson_errors=False)


plt.show()
#
# transport_times = []
# weights = []
# tube_lengths = []
# mean_trans_times = []
# mean_weights = []
# __tube_lenghts = []
#
# ax = plt.subplot()
# plt.figure()
# aux_ax = plt.subplot()
# for num in shots:
#
#     shot = Shot(num)
#     if shot.shotnum == 122:
#         shot.list.plot_erg_spectrum(200, 250)
#     sig, bg, bins = shot.list.get_time_dependence(218.5, signal_window_kev=1.5, nominal_values=True, bins=time_bins,
#                                                   normalization=600/shot.n_pulses, convolve=2/b_width)
#     # sig = convolve_gauss(sig, int(2/b_width))
#     bin_centers = (bins[:-1] + bins[1:])/2
#
#     def get_decay_cors(bs):
#         bl = bs[:-1]
#         br = bs[1:]
#         return 1.0/(0.5**(bl/40)-0.5**(br/40))
#
#     decay_corrs = get_decay_cors(bins)*1.0/shot.n_pulses*600
#     sig *= decay_corrs
#     arg_max = np.argmax(sig)
#     # sig = sig[: arg_max]
#     valid_points = np.where((1.0/np.sqrt(sig) < 0.05) & (sig < sig[arg_max]*0.7))
#     _weights = np.gradient(sig[valid_points])
#     _weights /= sum(_weights)
#     _weights = np.where(_weights >0, _weights, 1E-3)
#     __tube_lenghts.append(shot.tube_len)
#
#     times = bin_centers[valid_points]
#     mean_trans_times.append(np.average(times, weights=_weights))
#
#     mean_time = np.average(times, weights=_weights)
#     sigma = np.sqrt(np.average((mean_time-times)**2, weights=_weights))
#     mean_weights.append(1.0/sigma)
#     # transport_times.append(ufloat(mean_time, sigma))
#     # tube_lengths.append(shot.tube_len)
#
#     weights.extend(100*_weights)
#     transport_times.extend(times)
#     tube_lengths.extend([shot.tube_len]*len(times))
#     aux_ax.plot(times, _weights, label=shot.tube_len)
#
#     shot.list.plot_time_dependence(218.5, label=shot.__repr__(['tube_len']), ax=ax, bins=time_bins,
#                                    normalization=decay_corrs, convolve=6/b_width)
# aux_ax.legend()
# fig, ax = plt.subplots()
# m = LinearModel()
# fit = m.fit(data=mean_trans_times, x=__tube_lenghts, weights=np.array(mean_weights))
# # fit = LinearFit(__tube_lenghts, mean_trans_times)
# ax.scatter(tube_lengths, transport_times, s=weights)
# ax.scatter(__tube_lenghts, mean_trans_times, color='red', label='mean', marker='d', s=12)
# ax.legend()
#
# fit.plot_fit(ax=plt.subplots()[1])
#
#
# print(fit.fit_report())
# #
# print("velocity: ", 1.0/fit.params['slope'].value, 'm/s')
# #
# plt.show()


# print(shot2 == shot1)
# print(shot.tube_len)

