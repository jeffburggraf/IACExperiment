"""
For 0.5 total gas flow_pat with pattern 010010, the chamber clearing time is about 10 seconds
"""
import matplotlib.pyplot as plt
import numpy as np
from analysis import Shot, get_merged_time_dependence
# import numpy as np
from JSB_tools import mpl_hist, convolve_gauss
# import matplotlib.pyplot as plt
from JSB_tools.regression import LinearFit
from uncertainties import ufloat
from typing import Dict, List
import uncertainties.unumpy as unp
from JSB_tools.nuke_data_tools import Nuclide
from JSB_tools.nuke_data_tools.nuclide.fission_yields import FissionYields
from lmfit.models import LinearModel

nuclide_lines = {'Xe140': 806, "Xe139": 218.5, 'Nb99': 137.7, 'Sr94': 1427.7}
#  =============================================
shots_specification = {0.25: [9, 45, 46, 26, 27, 22, 23, 122],
                       0.5: [24, 25, 28, 29, 30, 31, 47, 48, 4, 5, 121, 120, 58, 59]}
shots = shots_specification[0.5]
# b_width = 5
flow = 0.25
nuclide = "Xe140"   #   805 for Xe-140
signal_window = 4
bg_window = 30
decay_corr = True  # needs work?
debug = False
debug_bins = np.linspace(0, 100, 5)
#  =============================================

# Shot(9).list.plotly(remove_background=True)
erg = nuclide_lines[nuclide]
shots_dict: Dict[float, List[Shot]] = {k: list(map(Shot, value)) for k, value in shots_specification.items()}
data = {}

nuclide = Nuclide(nuclide)
fiss_yields = FissionYields('U238', 'gamma')
print('Daughter yield: ', sum(fiss_yields.get_yield(nuclide.name)))

fiss_yields.plot_A()
for n_par in nuclide.get_decay_parents():  # print  fission yield data to consider parent decays
    p = f'Parent: {n_par} '
    try:
        p += f', fiss_yield: {np.mean(fiss_yields.yields[n_par.name])}'
        p += f' decay_modes: {n_par.decay_modes}'

    except KeyError:
        pass
    print(p)

hl = nuclide.half_life.n
max_time = 5*hl
b_width = hl/3


shots_dict[list(shots_dict.keys())[0]][-1].list.SPE.plot_erg_spectrum()

time_bins = np.arange(0, max_time + 1)

for shot in shots_dict[flow]:

    shot: Shot
    # time_dep, _, _ = shot.list.get_time_dependence(erg, debug_plot=False, bins=time_bins, signal_window_kev=window,
    #                                                bg_window_kev=13, nominal_values=False)
    time_dep, _, _ = shot.list.get_time_dependence(erg, debug_plot=False, bins=time_bins,
                                                   signal_window_kev=signal_window,
                                                   bg_window_kev=bg_window, nominal_values=True)

    if debug:
        if shot.shotnum == 23:
            print()
        shot.list.plot_time_dependence(erg, debug_plot=True, bins=debug_bins,
                                       signal_window_kev=signal_window,
                                       bg_window_kev=bg_window)
    time_dep = convolve_gauss(time_dep, b_width)
    time_dep /= shot.n_pulses
    # data[tube_length]: {'data': list[list[counts_vs_time]], 'shotnums': corresponding_shot_numbers}
    try:
        data[shot.tube_len]['data'].append(time_dep)
        data[shot.tube_len]['shotnums'].append(shot.shotnum)
    except KeyError:
        data[shot.tube_len] = {'data': [time_dep], 'shotnums': [shot.shotnum]}


x = []
y = []
yerr = []
fig, ax = plt.subplots()
for length, data_dict in data.items():
    data = data_dict['data']
    shotnums = data_dict['shotnums']
    sums = np.sum(data, axis=0)  # sum of counts for each time bin over all shots with same config
    stds = np.std(data, axis=0)  # STD of counts for each time

    sums /= max(sums)  # normalize to max of 1
    stds /= max(stds)

    # plot total time dependence of all configurations
    mpl_hist(time_bins, sums, ax=ax, label=length, poisson_errors=False, title='Mean counts vs time')

    select = np.where(np.array([max(q) for q in data]) > 0)  # max is a negative number: indicative of unusable statistics
    data = np.array(data)[select]  # Select shots
    shotnums = np.array(shotnums)[select]  # Select shot nums

    if decay_corr:  # correct for fact that x% of nuclei have decayed in each time bin
        _argmax = max([np.argmax(x) for x in data])  # Will be used to cut data off after max
        data_stnd = np.array([d[:_argmax+1] for d in data])
        corr = 1.0/np.array([0.5 ** (b1 / hl) - 0.5 ** (b2 / hl) for b1, b2 in
                             zip(time_bins[:_argmax+1], time_bins[1:_argmax+2])])  # decay correction
        data_stnd *= corr
        data_stnd = np.array([q/max(q) for q in data_stnd])  # standardize to max of data = 1
    else:
        _argmax = len(data[0])  # use 0 bs all lists in dat are same length
        data_stnd = np.array([q/max(q) for q in data])  # standardize to max of data = 1

    half_times = [time_bins[np.where(_c >= 0.5)][0] for _c in data_stnd]  # time in which half mx is reached
    x.append(length)  # For length vs median transport time
    y.append(np.mean(half_times))  # For length vs median transport time
    if len(half_times) > 1:
        yerr.append(np.std(half_times))
    else:
        yerr.append(0.1*np.mean(half_times))  # dont take standard dev of a single shot


    _ts = (time_bins[1:_argmax+2] + time_bins[:_argmax+1])/2
    plt.figure()
    for shotnum, d in zip(shotnums, data_stnd):  # plot individual shots
        # d /= max(d)
        plt.plot(_ts, d, label=shotnum)
        if decay_corr:
            plt.ylabel("Decay corrected counts")
        else:
            plt.ylabel("counts")

    plt.legend()
    plt.title(f'Tube length: {length}')

yerr = np.array(yerr)
plt.figure()
plt.errorbar(x, y, yerr, ls='None', marker='o')
plt.xlabel("Tube length [m]")
plt.ylabel(r"$\Delta$t to $\frac{1}{2}$ max count rate [s]")

model = LinearModel()
params = model.make_params()
weights = np.where(yerr != 0, yerr, 0.1)
weights = 1/weights
fit = model.fit(data=y, x=x, weights=weights,  params=params)
params = fit.params
eval_x = np.linspace(0, max(x), 100)
print(params)
plt.plot(eval_x, fit.eval(params=params, x=eval_x),
         label=f'lin fit, intercept =  {params["intercept"].value:.1f} +/- {params["intercept"].stderr:.1f}'
               f' slope =  {params["slope"].value:.1f} +/- {params["slope"].stderr:.1f}')
plt.legend()

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

