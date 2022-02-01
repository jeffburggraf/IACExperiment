import pickle

import matplotlib.pyplot as plt
from JSB_tools import Nuclide
from JSB_tools.spe_reader import SPEFile
from pathlib import Path
import re
from openpyxl import load_workbook
from cal_sources import CalSource
from scipy.signal import find_peaks
import uncertainties.unumpy as unp
import numpy as np
from JSB_tools.regression import LogPolyFit, PeakFit


def recal(spe: SPEFile, nuclide_name, intensity_thresh=0.03):
    nuclide = Nuclide.from_symbol(nuclide_name)
    peaks_ids, _ = find_peaks(unp.nominal_values(spe.counts), prominence=10*unp.std_devs(spe.counts) )
    peak_centers = spe.energies[peaks_ids]
    # ax = spe.plot_erg_spectrum(remove_baseline=True)
    _max = max(unp.nominal_values(spe.counts))
    g_lines = [g for g in nuclide.decay_gamma_lines if g.intensity >= intensity_thresh]
    true_ergs = [g.erg.n for g in g_lines]
    peak_channels = [peaks_ids[np.argmin(np.abs(erg-peak_centers))] for erg in true_ergs]
    meas_ergs = [peak_centers[np.argmin(np.abs(erg-peak_centers))] for erg in true_ergs]
    return np.polyfit(peak_channels, true_ergs, deg=2)[::-1]


# erg_cal = recal(SPEFile('/efficiencies/data/20211019-1723170-GEM20P4-70, 61-TP24472A.Spe'), 'Eu152')

data_path = Path(__file__).parent/'data'


wb = load_workbook(filename=data_path / 'Burggraf efficiency counts.xlsx')


def get_sources():
    sources = wb['Sources']

    headers = []
    for cell in sources['2']:
        headers.append(cell.value)

    cell_num = 3
    rows = []
    while True:
        row = {}
        for head, value in zip(headers, sources[str(cell_num)]):
            value = value.value
            row[head] = value
            # print(head, value)
        if row['EXPTID'] is None:
            break
        rows.append(row)
        cell_num += 1
    out = {}
    for row in rows:
        name = row['EXPTID'][2:]
        source_id = row['SAMPID']
        source = CalSource(name, source_id, row['Activity (Bq)'], row['Tzero(PST) (mm/dd/yy)'], unit='Bq')
        out[source_id] = source
    return out


def get_files():
    out = {}
    files = wb['Efficiency Counts']
    headers = [h.value for h in files['4']]
    headers[0] = 'stamp'
    for pairs in [(files['C'], files['F'], files['G']), (files['J'], files['M'], files['N'])]:
        for index, (src_id, dist, file_name) in enumerate(zip(*pairs)):
            src_id = src_id.value
            if index <2:
                continue
            if index > 20:break
            if src_id is None:
                continue
            f_path = data_path/file_name.value
            try:
                out[src_id][dist.value] = f_path
            except KeyError:
                out[src_id] = {dist.value: f_path}
    return out


sources = get_sources()
files = get_files()


#  ===============================================
threshold_intensity = 0.1
window = 2.85
debug_plot = True
#  ===============================================


x = {15.7: {}, 0.317: {}}
y = {15.7: {}, 0.317: {}}

effs = {15.7: {}, 0.317: {}}

ax_by_dist_and_src_name = {}

for src_id, dist_file_dict in get_files().items():
    for dist, fpath in dist_file_dict.items():
        source = CalSource.get_source(src_id)
        if source.name == 'CS137':
            nuclide = Nuclide.from_symbol("Ba137_m1")
        elif source.name == 'CD109':
            nuclide = Nuclide.from_symbol('Ag109_m1')
        else:
            nuclide = source.nuclide

        if source.name == 'BI207':
            window_kev = window*1.5
        else:
            window_kev = window

        gamma_lines = [g for g in nuclide.decay_gamma_lines if g.intensity > threshold_intensity or source.name == 'CD109']
        if source.name == 'TH228':
            gamma_lines = [Nuclide.from_symbol('Ra224').decay_gamma_lines[0]]
            window_kev = 1.9
        elif source.name == 'RA226':
            gamma_lines = [nuclide.decay_gamma_lines[0]]
        elif source.name == 'EU154':
            gamma_lines = gamma_lines[1:]

        if dist not in [0.317, 15.7]:
            scale = (dist/15.7)**2
        else:
            scale = 1

        if abs(dist - 15.70) < 1:
            dist = 15.70

        try:
            ax = ax_by_dist_and_src_name[source.name][dist]
        except KeyError:
            fig, axs = plt.subplots(1, 2)
            axs = axs.flatten()
            ax_by_dist_and_src_name[source.name] = {0.317: axs[0], 15.7: axs[1]}
            ax = ax_by_dist_and_src_name[source.name][dist]

        spe = SPEFile(fpath)
        n_decays = source.get_n_decays(spe.livetime, spe.system_start_time)
        ax.set_title(f'{source.name} @ {dist}')

        for g in gamma_lines:
            delta_window = 0 if g.erg < 750 else 2
            counts_meas, erg_bins = spe.get_counts(g.erg.n - window_kev/2 - delta_window,
                                                   g.erg.n + window_kev/2 + delta_window +0.4,
                                                   debug_plot=ax if debug_plot else False,
                                                   deadtime_corr=False, remove_baseline=True, make_rate=True,
                                                   return_bin_edges=True)
            counts_meas *= spe.livetime
            counts_meas = sum(counts_meas)

            true_counts = g.intensity*n_decays
            if source.name == 'CS137':
                true_counts *= 0.947
            eff = counts_meas/true_counts*scale
            ax.text(g.erg.n, ax.get_ylim()[-1]*0.5, f'eff = {eff:.1e}', rotation=60)
            n_label = source.nuclide.name
            try:
                effs[dist][n_label]['x'].append(g.erg.n)
                effs[dist][n_label]['y'].append(eff)
                # y[dist][n_label].append(eff)
            except KeyError:
                effs[dist][n_label] = {'x': [g.erg.n], 'y': [eff]}
                # y[dist][n_label] = [eff]

    # if debug_plot:
    #     plt.show()

fig, axs = plt.subplots(1, 2)
axs = axs.flatten()

xy_points = {0.317:  {'x': [], 'y': []}, 15.7: {'x': [], 'y': []}}
axes = {k: ax for k, ax in zip(effs.keys(), axs)}
for dist, labels_dict in effs.items():
    ax = axes[dist]
    ax.set_title(f"Distance  = {dist}")
    for nuclide_name, xy_dict in labels_dict.items():
        xs = xy_dict['x']
        ys = unp.nominal_values(xy_dict['y'])
        yerr = unp.std_devs(xy_dict['y'])
        ax.errorbar(xs, ys, yerr, label=nuclide_name, marker='o', ls='None')
        xy_points[dist]['x'].extend(xs)
        xy_points[dist]['y'].extend(xy_dict['y'])

    ax.legend()
    ax.set_xlabel("Energy")
    ax.set_ylabel("Efficiency")


_, fit_axs = plt.subplots(1, 2)

# fig, ax = plt.subplots()
_eval_x = np.linspace(0, 1000, 15)
_eval_y = []
for ax, (dist, xy_dist) in zip(fit_axs, xy_points.items()):
    ax.set_title("Log poly fits")
    ax.set_xlabel("Energy [Kev]")
    ax.set_ylabel("Efficiency")

    x = xy_dist['x']
    y = unp.nominal_values(xy_dist['y'])
    yerr = unp.std_devs(xy_dist['y'])

    fit_result = LogPolyFit(x, y, yerr, fix_coeffs=[0, 1, 2], order=4)
    f_path = '3_17mm' if dist == 0.317 else '15_5cm'
    with open(f"{f_path}.pickle", 'wb') as f:
        pickle.dump(fit_result.fit_result, f)

    fit_result.plot_fit(ax=ax)
    _eval_y.append(fit_result.fit_result.eval(x=_eval_x))
    print(f"Fit report for {dist}")
    print(fit_result.fit_result.fit_report())
    if dist == 15.7:
        x = list(sorted(x))
        print("Energy   Efficiency")
        for _x in x:
            print(f'{_x:<10.1f} {fit_result.eval_fit([_x])[0].n:.3e}')

print("Close/far ratio")
for i in range(len(_eval_x)):
    print(f'{_eval_x[i]} KeV: {_eval_y[0][i]/_eval_y[1][i]}')

plt.show()

