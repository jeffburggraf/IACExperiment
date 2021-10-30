import matplotlib.pyplot as plt

from JSB_tools.list_reader import MaestroListFile
import numpy as np
from lmfit.models import GaussianModel
from analysis import Shot
from scipy.signal import find_peaks
from uncertainties import unumpy as unp
from JSB_tools import Nuclide
from analysis import _get_maesto_list_shot_paths
from JSB_tools.list_reader import MaestroListFile
from JSB_tools.spe_reader import SPEFile

# s = Shot(121)
# s.llnl_spe.plot_erg_spectrum(remove_baseline=True)
# s.list.plotly(remove_baseline=True)

max_time = 230
shots = [[16, 25]]
peaks = [218.59, 397.4, 511, 602.35, 1427.7]
window_kev = 3
save_new_cal = True


def get_peaks(energies, counts):
    peak_ids, peak_infos = find_peaks(unp.nominal_values(counts), height=3*unp.std_devs(counts), width=2)
    fig, ax = plt.subplots()
    ax.plot(energies, unp.nominal_values(counts))
    ax.scatter(energies[peak_ids], peak_infos['peak_heights'], s=4, c='red')

    return ax, peak_ids, peak_infos


fig, all_ax = plt.subplots()


for shot_list in shots:
    counts = None
    erg_cals = None
    ergs = None
    erg_bin_widths = None
    channels = None
    list_files = []
    for shot_num in range(*shot_list):
        # shot = Shot(shot_num, load_erg_cal=False)
        try:
            list_file = MaestroListFile.from_pickle(_get_maesto_list_shot_paths()[shot_num]).SPE
        except KeyError:
            continue
        # list_file

        try:
            _counts = list_file.get_counts(remove_baseline=True)
            list_file.plot_erg_spectrum(ax=all_ax, leg_label=f"Shot{shot_num}", remove_baseline=True)
            list_files.append(list_file)
        except FileNotFoundError:
            continue

        if counts is None:
            ergs = list_file.energies
            counts = _counts
            erg_cals = list_file.erg_calibration
            erg_bin_widths = list_file.erg_bin_widths
            channels = list_file.channels
        else:
            counts += _counts
        if not all(c1 == c2 for c1, c2 in zip(list_file.erg_calibration, erg_cals)):
            print('Erg calib. change: ', list_file.erg_calibration, 'vs', erg_cals)
    ax, peaks_idx, _ = get_peaks(ergs, counts)
    calib_x, calib_y = [], peaks

    for peak in peaks:
        ax.plot([peak]*2, ax.get_ylim(), c='black')
        closest_idx = peaks_idx[np.argmin(np.abs(ergs[peaks_idx] - peak))]

        window = int(window_kev/(erg_bin_widths[closest_idx]))
        _slice = slice(closest_idx-window//2, closest_idx + window//2)
        fit_counts = counts[_slice]
        fit_channels = channels[_slice]
        model = GaussianModel()
        fit_y = unp.nominal_values(fit_counts)
        yerr = unp.std_devs(fit_counts)
        fit_weights = 1.0/yerr
        params = model.guess(data=fit_y, x=fit_channels)
        fit_result = model.fit(data=fit_y, x=fit_channels, weights=fit_weights, params=params)
        calib_x.append(params['center'])
        plt.figure()
        _ax = fit_result.plot_fit()
        _ax.set_title(peak)
    new_calib = list(reversed(np.polyfit(calib_x, calib_y, deg=2)))
    print(f'New cal: {new_calib}')
    print(f'Old cal {erg_cals}')
    new_ergs = np.sum([channels**i*c for i, c in enumerate(new_calib)], axis=0)
    ax.plot(new_ergs, unp.nominal_values(counts), label='New erg cal')
    ax.legend()
    if save_new_cal:
        for list_file in list_files:
            list_file.erg_calibration = new_calib
            list_file.pickle()
        # plt.plot(fit_channels, unp.nominal_values(fit_counts), label=f'{peak}')
plt.show()