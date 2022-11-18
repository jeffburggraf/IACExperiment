"""
Tests filter efficiency.
Results not super accurate due to poor knowledge of IAC efficiency and lack of rigor.
"""
import numpy as np
from matplotlib import pyplot as plt
from JSB_tools import mpl_hist
from pathlib import Path
from analysis import Shot
from JSB_tools import rebin, calc_background, shade_plot
from JSB_tools.nuke_data_tools import Nuclide
from JSB_tools.spe_reader import SPEFile
from mpant_reader import MPA
# =============================================
max_shots = 4
shot_kwargs = {'tube_len': None,  # 12.64, 9.44,  6.14, 4.16
               'foil_pos': 'upstream',  # 'center', 'upstream'
               'flow_pat': None,
               'he_flow': None,
               'ar_flow': None,
               'mylar': 0,
               'num_filters': 2,
               'flow_stop': 0,  # 0 for no stop
               'cold_filter': False
               }
e_window_width = 4
tests = [
    ('Tc101', 0),
    ('Tc104', 0),
    ('Cs140', 0),
    ('Sr94', 0)]

# =============================================

i = 0
y_iac = None
y_llnl = None
erg_bins = None


Shot.background_spe()
for shot in Shot.find_shots(**shot_kwargs):
    print(shot)
    i += 1

    if erg_bins is None:
        erg_bins = shot.list.erg_bins
        erg_bins = erg_bins[np.where((erg_bins > 50) & (erg_bins < 1600))]

    try:
        _y_iac = shot.iac_spe.counts/shot.iac_spe.effs
    except FileNotFoundError:
        continue

    _y_llnl = shot.list.get_erg_spectrum(eff_corr=True)
    _y_llnl /= shot.list.livetimes[-1]

    _y_iac = _y_iac / shot.iac_spe.livetime

    _y_llnl = rebin(shot.list.erg_bins, erg_bins, _y_llnl)
    _y_iac = rebin(shot.iac_spe.erg_bins, erg_bins, _y_iac)

    if y_llnl is None:
        y_llnl = _y_llnl
        y_iac = _y_iac
    else:
        y_llnl += _y_llnl
        y_iac += _y_iac

    if i >= max_shots:
        break


fig, (ax1, ax2) = plt.subplots(2, 1, sharex='all')

mpa_bg = MPA("/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/MCA/shot140.mpa")
mpa_bg.unpickle_eff("/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/MCA/cal_files/shot125_mpa.eff")

spe_bg = SPEFile("/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/shot140.Spe")
spe_bg.unpickle_eff("/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/cal_files/shot120.eff")

# spe_bg.pickle_eff()
# mpa_bg.pickle_eff()

bg_erg_llnl, bg_bins_llnl = spe_bg.get_counts(erg_min=60, erg_max=1600, make_rate=True, eff_corr=True,
                                              return_bin_edges=True)
bg_erg_iac, bg_bins_iac = mpa_bg.get_counts(erg_min=60, erg_max=1600, make_rate=True, eff_corr=True,
                                            return_bin_edges=True)

bg_erg_iac = rebin(bg_bins_iac, erg_bins, bg_erg_iac)
bg_erg_llnl = rebin(bg_bins_llnl, erg_bins, bg_erg_llnl)

mpl_hist(erg_bins, bg_erg_iac, label="IAC - BG", ax=ax2)
mpl_hist(erg_bins, bg_erg_llnl, label="LLNL - BG", ax=ax2)

y_iac = y_iac - calc_background(y_iac)
y_llnl = y_llnl - calc_background(y_llnl)

mpl_hist(erg_bins, y_iac, label="IAC", ax=ax1)
mpl_hist(erg_bins, y_llnl, ax=ax1, label="LLNL")

for (n_name, gi) in tests:
    erg = Nuclide(n_name).decay_gamma_lines[gi].erg.n
    emin = erg - e_window_width
    emax = erg + e_window_width
    i0 = np.searchsorted(erg_bins, emin)
    i1 = np.searchsorted(erg_bins, emax)
    tot_llnl = sum(y_llnl[i0: i1])
    tot_iac = sum(y_iac[i0: i1])
    print(f"{n_name}: {tot_llnl/(tot_iac + tot_llnl)}")
    ax1.axvline(erg)
    shade_plot(ax1, [emin, emax], label=n_name)


ax1.legend()

plt.show()
