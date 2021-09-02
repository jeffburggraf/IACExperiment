import warnings
from pathlib import Path
import re
from JSB_tools.spe_reader import SPEFile
import pylab as p
from JSB_tools.nuke_data_tools import GammaLine, Nuclide
from JSB_tools import mpl_hist, calc_background
from matplotlib import pyplot as plt
import numpy as np
from typing import Tuple
from JSB_tools.list_reader import MaestroListFile
from mpant_reader import MPA
from pickle_shots import get_maesto_list_shot_paths, get_mpant_mca_shot_paths
from shot_des import all_shots
import uncertainties.unumpy as unp
import matplotlib
from uncertainties import ufloat
# from lmfit import Parameters
shot_dirs = []
# matplotlib.use('TkAgg')
mca_paths = get_mpant_mca_shot_paths()
maestro_list_paths = get_maesto_list_shot_paths()


default_windows = {'llnl': {"Xe139-0": {'center': (216, 220.5), 'left': (200, 208), 'right': (223, 226)}}}


class Shot:
    def __init__(self, shot_num):
        self.mca: MPA
        self.shot_num = int(shot_num)
        try:
            self.mca = MPA(mca_paths[self.shot_num])
        except (FileNotFoundError, KeyError):
            warnings.warn(f"No MPA data for shot {shot_num}")
            self.mca: MPA = None
        try:
            self.maestro_list = MaestroListFile.from_pickle(maestro_list_paths[shot_num])
        except (FileNotFoundError, KeyError):
            raise

        self.maestro_zero_time = np.median(self.maestro_list.realtimes[np.where(self.maestro_list.sample_ready_state == 1)])
        self.maestro_list.times = self.maestro_list.times - self.maestro_zero_time
        self.maestro_list.__needs_updating__ = True

    def get_binned_times(self, time_bins, erg_min, erg_max, det='llnl') -> Tuple[np.ndarray, np.ndarray]:
        """
        Returns an array of the number of counts in each time bin subject to an energy cut.
        Args:
            time_bins:
            erg_min:
            erg_max:
            det:

        Returns: Bin values,  bin centers

        """
        time_bins = np.array(time_bins)
        if det.lower() == 'llnl':
            times = self.maestro_list.get_times_in_range(erg_min=erg_min, erg_max=erg_max)
            values, _ = np.histogram(times, bins=time_bins)
            b_centers = (time_bins[1:] + time_bins[:-1])/2
            return b_centers, unp.uarray(values, np.sqrt(values))
        else:
            raise NotImplementedError

    @staticmethod
    def get_rise_an_fall_times(_times, _weights):
        rise_time_i = np.argmax(_weights > max(_weights) / 2)
        rise_time = _times[rise_time_i]
        kernal = np.ones(4) / 4.
        _times = np.convolve(_times, kernal, mode='same')
        max_time_i = np.argmax(_weights)
        _weights = _weights[max_time_i:]
        _times = _times[max_time_i:] - _times[max_time_i]
        half_time_i = np.argmin(np.abs(_weights - _weights[0] / 2))
        # plt.plot(_times, _weights)

        return rise_time, _times[half_time_i]

    def plot_time_dependence(self, nuclide='Xe139', gamma_index=0,
                             det='llnl', time_bins=np.arange(0, 40*6, 2), subfig=False, verbose=False):
        window_entry = default_windows[det.lower()][f'{nuclide}-{gamma_index}']
        window_left = window_entry['left']
        window_right = window_entry['right']
        window_signal = window_entry['center']
        _, bg_left = self.get_binned_times(time_bins, *window_left, det=det)
        _, bg_right = self.get_binned_times(time_bins, *window_right, det=det)
        _, sig = self.get_binned_times(time_bins, *window_signal, det=det, )
        bg_left *= (window_signal[1]-window_signal[0])/(window_left[1]-window_left[0])
        bg_right *= (window_signal[1]-window_signal[0])/(window_right[1]-window_right[0])
        bg = (bg_left + bg_right)/2
        sig -= bg
        primary_trans_time, hl = self.get_rise_an_fall_times((time_bins[1:]+time_bins[:-1])/2, sig)
        if subfig:
            fig = plt.figure(constrained_layout=True)
            subfigs = fig.subfigures(nrows=2, ncols=1)
            ax1, ax2 = subfigs[0].subplots(nrows=1, ncols=2)
            ax3, ax4 = subfigs[1].subplots(nrows=1, ncols=2)
            subfigs[0].suptitle("Detector 1")
            subfigs[1].suptitle("Detector 2")
        else:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        fig.suptitle(f'shot {self.shot_num}')
        ax3.text(0.15, 0.15, 'NO DATA (YET)', transform=ax3.transAxes, fontsize=20, rotation=30)
        mpl_hist(time_bins, sig, ax=ax1, label='signal')
        ax1.text(0.6, 0.6, f'rise time:\n{primary_trans_time:.1f}')
        ax1.legend()
        llnl_tot_counts = np.sum(sig)
        ax2.text(0.7, 0.5, f'Tot events\n{llnl_tot_counts:.1e}', transform=ax2.transAxes)

        gamma_line = Nuclide.from_symbol(nuclide).decay_gamma_lines[gamma_index]
        erg = gamma_line.erg.n
        half_xrange = max(map(abs, [window_left[0]-erg, window_right[-1]-erg])) + 2

        self.maestro_list.plot_erg_spectrum(erg-half_xrange, erg+half_xrange, time_bins[0], time_bins[-1],
                                            ax=ax2)
        self.mca.plot_spectrum(ax=ax4, erg_min=erg-half_xrange, erg_max=erg+half_xrange)
        _lims = ax2.get_ylim()
        ax2.fill_between(window_right, [_lims[0]]*2, [_lims[1]]*2, alpha=0.3, color='red', label='bg window')
        ax2.fill_between(window_left, [_lims[0]]*2, [_lims[1]]*2, alpha=0.3, color='red')
        ax2.fill_between(window_signal, [_lims[0]]*2, [_lims[1]]*2, alpha=0.3, color='blue', label='signal window')
        ax2.legend()

        plt.subplots_adjust(hspace=0.37)
        iac_counts = self.mca.counts - calc_background(self.mca.get_counts(nominal=False))
        iac_tot_counts = np.sum(iac_counts[np.where((self.mca.energies >= window_signal[0]) &
                                                    (self.mca.energies <= window_signal[1]))])

        ax4.text(0.7, 0.5, f'Tot events\n{iac_tot_counts:.1e}', transform=ax4.transAxes)

        ax1.set_xlabel("Time since shot [s]")
        ax1.set_ylabel("counts")
        ax3.set_xlabel("Time since shot [s]")
        ax3.set_ylabel("counts")
        return ax1


if __name__ == '__main__':
    # s = Shot(134)
    # s = Shot(139)
    s = Shot(134)
    spe = SPEFile('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/tuesday/BG.Spe')

    sig, _ = np.histogram(s.maestro_list.energies, s.maestro_list.erg_bins)
    sig = sig / s.maestro_list.livetimes[-1]
    sig = unp.uarray(sig, np.sqrt(sig))
    sig = sig - calc_background(sig)
    bg = spe.counts
    bg = bg / spe.livetime
    bg = bg - calc_background(bg)
    ax, _ = mpl_hist(s.maestro_list.erg_bins, sig, label="Signal")
    mpl_hist(spe.erg_bins, 5*bg, ax=ax, label='Background')
    ax.legend()
    ax.set_xscale('log')
    ax.set_xlabel('Energy')
    ax.set_ylabel('Counts/sec')
    # ax.set_yscale('log')

    # s.plot_time_dependence()

    # s.mca.plot_spectrum()

    # s = Shot(30)
    # s.mca.plot_spectrum()
    # MPA.()
    # s = Shot(134)
    # s.plot_time_dependence()
    # t_bins = np.arange(0, 200, 3)
    # x, y = s.get_binned_times(t_bins, 217, 220)
    # x1, y1 = s.get_binned_times(t_bins, 214, 217)
    #
    # plt.plot(x, y)
    # plt.plot(x1, y1, label='off-peak')
    # plt.legend()
    #
    # s.maestro_list.plot_erg_spectrum()

    plt.show()