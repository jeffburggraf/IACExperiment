"""
Reads data taken at the  LLNL counting facility and pickles dictionary of calculated efficiency of
each relevant gamma line (xy_points.pickle).
"""
import pickle
import warnings
import matplotlib.pyplot as plt
from JSB_tools import TabPlot
from JSB_tools.nuke_data_tools import Nuclide
from JSB_tools.spe_reader import SPEFile
from pathlib import Path
from openpyxl import load_workbook
from cal_sources import CalSource
import uncertainties.unumpy as unp
import numpy as np
from JSB_tools.regression import LogPolyFit


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


def get_gamma_lines_and_window(source):
    if source.name == 'CS137':
        nuclide = Nuclide("Ba137_m1")
    elif source.name == 'CD109':
        nuclide = Nuclide('Ag109_m1')
    else:
        nuclide = source.nuclide

    if source.name == 'BI207':
        window_kev = window * 1.5
    else:
        window_kev = window

    gamma_lines = [g for g in nuclide.decay_gamma_lines if g.intensity > threshold_intensity or source.name == 'CD109']
    if source.name == 'TH228':
        gamma_lines = [Nuclide('Ra224').decay_gamma_lines[0]]
        window_kev = 1.9
    elif source.name == 'RA226':
        gamma_lines = [nuclide.decay_gamma_lines[0]]
    elif source.name == 'EU154':
        gamma_lines = gamma_lines[1:]

    gamma_lines = [g for g in gamma_lines if (g.erg.n > erg_min and (g.parent_nuclide.name, round(g.erg.n, 2)) not in exception_gamma_lines)]

    return gamma_lines, window_kev


sources = get_sources()  # mapping of src_if to CalSource objs
files = get_files()  # files[src_id][dist] gives path to corresponding file


#  ===============================================
threshold_intensity = 0.02
window = 2.85
erg_min = 55
debug_plot = True
only_close = False
pickle_xypoints = True
exception_gamma_lines = [('Ba133', 79.61)]  # Excluded gamma lines. Round to two decimals.
#  ===============================================

# exception_gamma_lines = [Nuclide(n).get_gamma_nearest(erg) for n, erg in exception_gamma_lines]
print("Excluded gamma lines:")
for g in exception_gamma_lines:
    print(f"\t{g}")

effs = {0.317: {}, 15.7: {}}
xy_points = {0.317:  {'x': [], 'y': [], 'nuclides': []}, 15.7: {'x': [], 'y': [], 'nuclides': []}}

tab_fars = []

tab_close = TabPlot()

if not only_close:
    tab_fars.append(TabPlot())


for src_id, dist_file_dict in get_files().items():
    for dist, fpath in dist_file_dict.items():
        if only_close and dist != 0.317:
            continue

        source = CalSource.get_source(src_id)

        gamma_lines, window_kev = get_gamma_lines_and_window(source)

        if dist not in [0.317, 15.7]:
            scale = (dist/15.7)**2
        else:
            scale = 1

        if abs(dist - 15.70) < 1:
            dist = 15.70

        spe = SPEFile(fpath)
        assert spe.eff_path is None
        n_decays = source.get_n_decays(spe.livetime, spe.system_start_time)

        for g in gamma_lines:
            print(g)
            tab = None
            if dist == 0.317:
                tab = tab_close
            else:
                if not only_close:
                    if tab_fars[-1].max_buttons_reached:
                        tab_fars.append(TabPlot())
                    tab = tab_fars[-1]

            ax_spectra = None
            if tab is not None:
                ax_spectra, _ = tab.new_ax(f"{source.nuclide.name} {g.erg.n:.1f}", 1, 2)

            delta_window = 0 if g.erg < 750 else 2
            counts_meas, erg_bins = spe.get_counts(g.erg.n - window_kev/2 - delta_window,
                                                   g.erg.n + window_kev/2 + delta_window +0.4,
                                                   debug_plot=ax_spectra if debug_plot else False,
                                                   deadtime_corr=False, remove_baseline=True,
                                                   return_bin_edges=True)
            counts_meas = sum(counts_meas)

            true_counts = g.intensity*n_decays
            if source.name == 'CS137':
                true_counts *= 0.947  # Feeding ratio

            eff = counts_meas/true_counts*scale

            n_label = source.nuclide.name
            xy_points[dist]['x'].append(g.erg.n)
            xy_points[dist]['y'].append(eff)
            xy_points[dist]['nuclides'].append(n_label)

            try:
                effs[dist][n_label]['x'].append(g.erg.n)
                effs[dist][n_label]['y'].append(eff)
                # y[dist][n_label].append(eff)
            except KeyError:
                effs[dist][n_label] = {'x': [g.erg.n], 'y': [eff]}
                # y[dist][n_label] = [eff]

for d in xy_points.values():
    arg_srt = np.argsort(d['x'])
    d['x'] = np.array(d['x'])[arg_srt]
    d['y'] = np.array(d['y'])[arg_srt]
    d['nuclides'] = np.array(d['nuclides'])[arg_srt]


for tabs, dist in [([tab_close], 0.317), (tab_fars, 15.7)]:
    x = xy_points[dist]['x']
    y = xy_points[dist]['y']
    err = unp.std_devs(y)
    y = unp.nominal_values(y)

    model = LogPolyFit(x, y, err, fix_coeffs='all', order=4)
    print(f"{dist} cm model coeffs:\n\t{model.coefs}")

    for tab in tabs:
        for _, ax_persistent in tab.plt_axs:

            for label, xy in effs[dist].items():
                xs = xy['x']
                ys = xy['y']
                ax_persistent.errorbar(xs, unp.nominal_values(ys), unp.std_devs(ys), label=label, ls='None', marker='.',
                                       capsize=5, markersize=8)
                # ax_persistent.plot()
            ax_persistent.legend()

if pickle_xypoints:
    if not only_close:
        with open("xy_points.pickle", 'wb') as f:
            pickle.dump(xy_points, f)
    else:
        warnings.warn("Did not pickle xy_points bc only_close is True")

plt.show()

    # if debug_plot:
    #     plt.show()

# fig, axs = plt.subplots(1, 2)
# axs = axs.flatten()
#
#
#
# axes = {k: ax for k, ax in zip(effs.keys(), axs)}
# for dist, labels_dict in effs.items():
#     ax = axes[dist]
#     ax.set_title(f"Distance  = {dist}")
#     for nuclide_name, xy_dict in labels_dict.items():
#         xs = xy_dict['x']
#         ys = unp.nominal_values(xy_dict['y'])
#         yerr = unp.std_devs(xy_dict['y'])
#         ax.errorbar(xs, ys, yerr, label=nuclide_name, marker='o', ls='None')
#         xy_points[dist]['x'].extend(xs)
#         xy_points[dist]['y'].extend(xy_dict['y'])
#
#     ax.legend()
#     ax.set_xlabel("Energy")
#     ax.set_ylabel("Efficiency")
#
# for k, xy_dict in xy_points.items():
#     x, y, = map(np.array,  [xy_dict['x'], xy_dict['y']])
#     arg_srt = np.argsort(x)
#     xy_dict['x'] = x[arg_srt]
#     xy_dict['y'] = y[arg_srt]
#
# with open("xy_points.pickle", 'wb') as f:
#     pickle.dump(xy_points, f)
#
# _, fit_axs = plt.subplots(1, 2)
#
# # fig, ax = plt.subplots()
# _eval_x = np.linspace(0, 1000, 15)
# _eval_y = []
# for ax, (dist, xy_dist) in zip(fit_axs, xy_points.items()):
#     ax.set_title("Log poly fits")
#     ax.set_xlabel("Energy [Kev]")
#     ax.set_ylabel("Efficiency")
#
#     x = xy_dist['x']
#     y = unp.nominal_values(xy_dist['y'])
#     yerr = unp.std_devs(xy_dist['y'])
#
#     fit_result = LogPolyFit(x, y, yerr, fix_coeffs='all', order=5)
#     f_path = '3_17mm' if dist == 0.317 else '15_5cm'
#     with open(f"{f_path}.pickle", 'wb') as f:
#         pickle.dump(fit_result.fit_result, f)
#
#     fit_result.plot_fit(ax=ax)
#     _eval_y.append(fit_result.fit_result.eval(x=_eval_x))
#     print(f"Fit report for {dist}")
#     print(fit_result.fit_result.fit_report())
#     if dist == 15.7:
#         x = list(sorted(x))
#         print("Energy   Efficiency")
#         for _x in x:
#             print(f'{_x:<10.1f} {fit_result.eval_fit([_x])[0].n:.3e}')
#
# print("Close/far ratio")
# for i in range(len(_eval_x)):
#     print(f'{_eval_x[i]} KeV: {_eval_y[0][i]/_eval_y[1][i]}')
#
# plt.show()

