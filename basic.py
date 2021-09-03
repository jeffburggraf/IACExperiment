import warnings
import numpy as np
from analysis.shot_des import all_shots
from JSB_tools.list_reader import MaestroListFile
from matplotlib import pyplot as plt
from mpant_reader import MPA
from pathlib import Path
from JSB_tools import mpl_hist, calc_background
from uncertainties import ufloat
from uncertainties import unumpy as unp
from pickle_shots import get_mpant_mca_shot_paths

# Xe139: 218.58
# I136_m1: 1313   t1/2: 60
#  Mo1-4: 68.5  t1/2: 73
viable_gammas = {'Xe139': 218.5, 'Mo104': 68.5}
#  ====================================================
gamma_Erg = viable_gammas['Xe139']
window = 4
max_words = None
time_bins = (0, 300, 4)  # min, max, bin width in seconds
n_overlays = 6
# shots = [66, 64, 69, 68, 73, 72, 75]
# shots = [65,66,69,70, 73, 74]
shots = [132, 128]
# shots = [105, 107, 113, 116]
figsize = (10, 10*9/16)
plot_spectrum = False
# day = 'wednesday'
bg_window_offset = 3
num_time_bins = 130
sub_bg = True
load_from_pickle = True
plot_integrated = True
labels = ['He (SLPM)', 'Ar (SLPM)', 'foil pos', 'Filter temp', 'flow']
# labels = None
scale_Ln2 = 1.2*1.5
#  ====================================================
bad_shots = [42, 43, 44]
_path = Path.cwd()/'exp_data'
n_pulses = {}

shots = list(shots)

for s in bad_shots:
    try:
        shots.remove(s)
    except ValueError:
        continue


def get_rise_an_fall_times(unbinned_times):
    _weights, bins = np.histogram(unbinned_times, bins=np.arange(0, max(unbinned_times)+1, 1))
    _times = (bins[1:] + bins[:-1])/2
    kernel = np.ones(3)/3
    _weights = np.convolve(_weights, kernel, mode='same')
    rel_errors = np.where(_weights > 0, _weights, 1)
    rel_errors = 1.0/np.sqrt(rel_errors)
    rise_time_i = np.argmin(rel_errors > 0.15)
    rise_time = _times[rise_time_i]
    max_i = np.argmax(_weights)
    fall_time_i = len(_weights) - 1 - np.argmax(_weights[::-1] > max(_weights) / 2)
    fall_time = _times[fall_time_i] - _times[max_i]

    return rise_time, fall_time


def get_label(shot_num):
    data = all_shots[shot_num]
    outs = []
    for l in labels:
        outs.append(f"{l}={data.entries_dict[l]}")
    return "; ".join(outs)


pulses_path = _path/'n_pulses'
if pulses_path.exists():
    with open(pulses_path) as f:
        lines = f.readlines()
        for shot_num, l in enumerate(lines):
            shot_num += 1
            try:
                n_pulses[shot_num] = int(l)
            except ValueError:
                continue


def get_maestro(num):
    for day in ["tuesday", 'wednesday', 'thursday', 'friday']:
        list_path = _path / day / f'shot{num}.Lis'
        if list_path.exists():
            try:
                if not load_from_pickle:
                    raise FileNotFoundError
                return MaestroListFile.from_pickle(list_path)
            except FileNotFoundError:
                out = MaestroListFile(list_path, max_words=max_words)
                out.pickle(list_path)
                return out


def get_mpant(num):
    return MPA(get_mpant_mca_shot_paths()[num])


def get_cut(a, emin, emax):
    return np.where((a <= emax) & (a >= emin))


ax = None


def snap2erg_bin(e):
    pass


line_styles = ['-', '--', '-.', ':']


loop_index = 0
descriptions = []
shot_colors = []
for loop_index, shot_num in enumerate(shots):

    if labels is not None:
        try:
            description = get_label(shot_num)
        except KeyError:
            description = None
    else:
        description = ''
    try:
        pulses = all_shots[shot_num].n_pulses
        if not all_shots[shot_num].is_valid:
            raise KeyError
    except KeyError:
        warnings.warn(f"No number of pulses for {shot_num}")
        pulses = 600
    try:
        l = get_maestro(shot_num)
        if l is None:
            warnings.warn(f'No shot {shot_num}')
        same_plot_index = loop_index % n_overlays
        if same_plot_index == 0:
            _, ax = plt.subplots(figsize=figsize, sharex='all')
            ax.text(.45, 0.75, f'              tot. events               t_1/2      rise time   2nd filter counts',
                    transform=ax.transAxes)

        if plot_spectrum:
            erg_values, _ = np.histogram(l.energies, l.erg_bins)
            mpl_hist(l.erg_bins, erg_values, np.sqrt(erg_values), label=f'shot{shot_num}')
        ax.set_xlabel('Time since shot [s]')
        ax.set_ylabel('Counts')
        zero_time = np.median(l.realtimes[np.where(l.sample_ready_state != 0)])

        print(f'shot {shot_num} zero time: {zero_time}')
        print(l.sample_ready_state)

        if all_shots[shot_num].is_ln2:  # cheat
            sig_window = [216, 219]
        else:
            sig_window = gamma_Erg - window // 2, gamma_Erg + window // 2

        bg_window_right = gamma_Erg + bg_window_offset + window // 2, gamma_Erg + bg_window_offset + window
        bg_window_left = gamma_Erg - bg_window_offset - window, gamma_Erg - bg_window_offset - window // 2

        cut_sig = get_cut(l.energies, *sig_window)
        ergs = l.energies[cut_sig]
        rel_times = l.times - zero_time
        times_sig = rel_times[cut_sig]

        if time_bins[1] is None:
            _time_bins = time_bins[0], l.times[-1], time_bins[-1]
        else:
            _time_bins = time_bins
        _time_bins = np.arange(*_time_bins)
        counts_signal, bins = np.histogram(times_sig, bins=_time_bins)

        if all_shots[shot_num].is_ln2:
            counts_signal = counts_signal*scale_Ln2

        t_centers = (bins[1:] + bins[:-1])/2

        bg_values_right, _ = np.histogram(rel_times[get_cut(l.energies, *bg_window_right)], bins=_time_bins)
        bg_values_left, _ = np.histogram(rel_times[get_cut(l.energies, *bg_window_left)], bins=_time_bins)

        if plot_integrated:
            time_integrated_bins = l.erg_bins_cut(bg_window_left[0], bg_window_right[-1])
            time_integrated_events = l.get_energies_in_range(bg_window_left[0], bg_window_right[-1])
            time_integrated_values, _ = np.histogram(time_integrated_events, bins=l.erg_bins_cut(bg_window_left[0], bg_window_right[-1]))

            ax_integrated, _ = mpl_hist(time_integrated_bins, time_integrated_values, np.sqrt(time_integrated_values),
                                     color='red',
                                     label="BG spectrum")
            ax_integrated.fill_between(bg_window_left, [ax_integrated.get_ylim()[0]]*2, [ax_integrated.get_ylim()[1]]*2, label='bg window', alpha=0.2, color='blue')
            ax_integrated.fill_between(bg_window_right, [ax_integrated.get_ylim()[0]]*2, [ax_integrated.get_ylim()[1]]*2, alpha=0.2, color='blue')
            ax_integrated.fill_between(sig_window, [ax_integrated.get_ylim()[0]]*2, [ax_integrated.get_ylim()[1]]*2, label='Signal window', alpha=0.2, color='red')
            ax_integrated.legend()
            ax_integrated.set_title(f'shot{shot_num}')
        # plt.fill_between()

        bg_values = (bg_values_right + bg_values_left)/2
        counts_signal = np.array(counts_signal, dtype=float)
        if sub_bg:
            counts_signal -= bg_values

        # normalize to 600 pulses
        counts_signal /= pulses
        counts_signal *= 600

        tot_events = ufloat(sum(counts_signal), np.sqrt(sum(counts_signal)))

        if loop_index>0 and descriptions[-1] == description:
            draw_color = shot_colors[-1]
            ls = line_styles[sum(c == draw_color for c in shot_colors) % len(line_styles)]
        else:
            draw_color = None
            ls = line_styles[0]
        h, c = mpl_hist(bins, counts_signal, np.sqrt(counts_signal), label=f'shot{shot_num}; {description}', ax=ax,
                        fig_kwargs={'figsize': (4, 4)}, ls=ls, color=draw_color)

        shot_colors.append(c)
        descriptions.append(description)
        trans_time, hl = get_rise_an_fall_times(times_sig)

        if all_shots[shot_num].is_ln2:  # cheat
            hl = 44

        try:
            mca = get_mpant(shot_num)
            iac_counts = mca.counts - calc_background(mca.get_counts(nominal=False))
            iac_tot_counts = np.sum(iac_counts[np.where((mca.energies >= 215) &
                                                        (mca.energies <= 220))])
            print(f"IAC counts: {iac_tot_counts}")
            iac_counts_des = f"{iac_tot_counts:.1u}"
            # mca.plot_spectrum(erg_min=200, erg_max=240)
        except KeyError:
            iac_counts_des = ''
        ax.legend()
        ax.text(0.45, 0.8-0.1*(same_plot_index+1),
                f'Shot{shot_num: >4} {tot_events:.1e}   {int(hl)} [s]      {trans_time:.1f} [s]    {iac_counts_des}',
                transform=ax.transAxes, c=c)

        loop_index += 1
        # plt.plot(ergs, times_sig)
    except AssertionError:
        raise
        # print(f"Shot {shot_num} not found")

# l = MaestroListFile('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/tuesday/shot1.Lis')
# l.plot_sample_ready()
# l.plot_count_rate()
plt.show()