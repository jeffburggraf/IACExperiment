import warnings
import numpy as np
from exp_data.shot_des import all_shots
from JSB_tools.list_reader import MaestroListFile
from matplotlib import pyplot as plt
from pathlib import Path
import numpy as np
from JSB_tools import mpl_hist
from uncertainties import ufloat
from uncertainties import unumpy as unp

# Xe139: 218.58
# I136_m1: 1313   t1/2: 60
#  Mo1-4: 68.5  t1/2: 73
viable_gammas = {'Xe139': 218.5, 'Mo104': 68.5}
#  ====================================================
gamma_Erg = viable_gammas['Xe139']
window = 2
max_words = None
t_max = None
n_overlays = 2
# shots = [66, 64, 69, 68, 73, 72, 75]
shots = [125, 129, 130, 131, 132, 134]
figsize = (7, 5) #(10, 10*9/16)
plot_spectrum = False
# day = 'wednesday'
bg_window_offset = 3
num_time_bins = 130
sub_bg = True
load_from_pickle = True
plot_integrated = False
labels = ['He (SLPM)', 'Ar (SLPM)', 'Filter temp', 'Mylar (um)', 'flow']
# labels = None
#  ====================================================
_path = Path.cwd()/'exp_data'
n_pulses = {}


def get_rise_an_fall_times(_times, _weights):
    rise_time_i = np.argmax(_weights>max(_weights)/2)
    rise_time = _times[rise_time_i]
    kernal = np.ones(4)/4.
    _times = np.convolve(_times, kernal, mode='same')
    max_time_i = np.argmax(_weights)
    _weights = _weights[max_time_i:]
    _times = _times[max_time_i:] - _times[max_time_i]
    half_time_i = np.argmin(np.abs(_weights - _weights[0]/2))
    # plt.plot(_times, _weights)

    return rise_time, _times[half_time_i]


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


def get_cut(a, emin, emax):
    return np.where((a <= emax) & (a >= emin))


ax = None


def snap2erg_bin(e):
    pass


loop_index = 0
for shot_num in shots:

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
            ax.text(.45, 0.9, f'              tot. events                t_1/2    trans. time', transform=ax.transAxes)

        if plot_spectrum:
            erg_values, _ = np.histogram(l.energies, l.erg_bins)
            mpl_hist(l.erg_bins, erg_values, np.sqrt(erg_values), label=f'shot{shot_num}')

        zero_time = np.median(l.realtimes[np.where(l.sample_ready_state != 0)])
        print(f'shot {shot_num} zero time: {zero_time}')
        print(l.sample_ready_state)

        bg_window_right = gamma_Erg + bg_window_offset + window // 2, gamma_Erg + bg_window_offset + window
        bg_window_left = gamma_Erg - bg_window_offset - window, gamma_Erg - bg_window_offset - window // 2
        sig_window = gamma_Erg-window//2, gamma_Erg+window//2

        cut_sig = get_cut(l.energies, *sig_window)
        ergs = l.energies[cut_sig]
        rel_times = l.times - zero_time
        times_sig = rel_times[cut_sig]
        if t_max is None:
            t_max = l.times[-1]
        counts_signal, bins = np.histogram(times_sig, bins=np.linspace(0, t_max, num_time_bins))
        t_centers = (bins[1:] + bins[:-1])/2

        bg_values_right, _ = np.histogram(rel_times[get_cut(l.energies, *bg_window_right)], bins=np.linspace(0, t_max, num_time_bins))
        bg_values_left, _ = np.histogram(rel_times[get_cut(l.energies, *bg_window_left)], bins=np.linspace(0, t_max, num_time_bins))

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
        counts_signal = np.array(counts_signal, dtype=np.float)
        if sub_bg:
            counts_signal -= bg_values

        # normalize to 600 pulses
        counts_signal /= pulses
        counts_signal *= 600

        tot_events = ufloat(sum(counts_signal), np.sqrt(sum(counts_signal)))

        h, c = mpl_hist(bins, counts_signal, np.sqrt(counts_signal), label=f'shot{shot_num}; {description}', ax=ax, fig_kwargs={'figsize':(4, 4)})
        trans_time, hl = get_rise_an_fall_times(t_centers, counts_signal)

        ax.legend()
        ax.text(0.45, 0.9-0.1*(same_plot_index+1),
                f'Shot {shot_num} {tot_events}   {int(hl)} [s]    {trans_time:.1f} [s]',
                transform=ax.transAxes, c=c)

        loop_index += 1
        # plt.plot(ergs, times_sig)
    except AssertionError:
        print(f"Shot {shot_num} not found")

# l = MaestroListFile('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/tuesday/shot1.Lis')
# l.plot_sample_ready()
# l.plot_count_rate()
plt.show()