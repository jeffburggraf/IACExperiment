from pathlib import Path
import numpy as np
import re
from matplotlib import pyplot as plt
from datetime import datetime
from scipy.stats import norm
from functools import cached_property
from JSB_tools import mpl_hist
from uncertainties.core import AffineScalarFunc
from uncertainties import unumpy as unp
from typing import List


class MPA:
    def __init__(self, path):
        with open(path) as f:
            lines = f.readlines()
        header = lines[:lines.index("[DATA0,8192 ]\n")]
        adc_headers = []
        try:
            for i in range(6):
                start_i = header.index(f"[ADC{i+1}]\n")
                end_s = f"[ADC{i+2}]\n"
                if end_s not in header:
                    adc_headers.append(header[start_i+1:])
                else:
                    end_i = header.index(end_s)
                    adc_headers.append(header[start_i+1: end_i])
        except ValueError:
            pass
        adc_headers_dict = {}
        for index, head in enumerate(adc_headers):
            index += 1
            d = {}
            for s in head:
                try:
                    k, v = s.split('=')
                    v = v.rstrip().lstrip()
                    try:
                        v = float(v)
                    except ValueError:
                        pass
                    d[k] = v

                except ValueError:
                    pass
            if d['active']:
                adc_headers_dict[index] = d

        counts = [[] for _ in range(len(adc_headers_dict))]

        index = 0
        for line in lines[len(header)+1:]:
            if re.match("\[DATA[0-9],", line):
                index += 1
            else:
                counts[index].append(int(line))

        self._counts = [unp.uarray(c, np.sqrt(c)) for c in counts]

        channels = [np.arange(0, adc_dict['range']) + 0.5 for adc_dict in adc_headers_dict.values()]

        channels_bins = [np.arange(0, adc_dict['range']+1) for adc_dict in adc_headers_dict.values()]

        self._erg_coeffs = np.array([[adc_dict['caloff'], adc_dict['calfact']] for adc_dict in adc_headers_dict.values()])

        self._energies = [np.sum([c*channels[adc_i]**pow for pow, c in enumerate(self._erg_coeffs[adc_i])], axis=0)
                          for adc_i in range(len(channels))]
        self._erg_bins = [np.sum([c*channels_bins[adc_i]**pow for pow, c in enumerate(self._erg_coeffs[adc_i])], axis=0)
                          for adc_i in range(len(channels))]

        self._live_times = [adc_dict['realtime'] for adc_dict in adc_headers_dict.values()]
        self._real_times = [adc_dict['runtime'] for adc_dict in adc_headers_dict.values()]

        for line in header:
            if m := re.match("REPORT-FILE from (.+) written", line):
                s = m.groups()[0].strip()
                self.system_start_time = datetime.strptime(s, "%m/%d/%Y %H:%M:%S")
                break

    def channel2erg(self, ch, adc=1):
        # todo: Is this right? is there a "zero" channel? Compare to MPANT.
        return self.get_energies(adc)[int(ch)]

    @property
    def livetime(self):
        return self.get_livetime()

    def get_livetime(self, adc=1):
        return self._live_times[adc-1]

    @property
    def realtime(self):
        return self.get_realtime()

    def get_realtime(self, adc=1):
        assert adc > 0, "ADC starts at 1, not 0."
        return self._real_times[adc-1]

    @property
    def counts(self, adc=1):
        return self.get_counts(adc)

    def get_counts(self, adc=1):
        assert adc > 0, "ADC starts at 1, not 0."
        return self._counts[adc-1]

    @property
    def energies(self, adc=1):
        return self.get_energies(adc)

    def get_energies(self, adc=1):
        assert adc > 0, "ADC starts at 1, not 0."
        return self._energies[adc - 1]

    @property
    def erg_bins(self):
        return self.get_erg_bins(adc=1)

    def get_erg_bins(self, adc=1):
        assert adc > 0 , "ADC starts at 1, not 0."
        return self._erg_bins[adc-1]

    def erg_bin_index(self, erg):
        if hasattr(erg, '__iter__'):
            return np.array([self.erg_bin_index(e) for e in erg])
        if isinstance(erg, AffineScalarFunc):
            erg = erg.n
        return np.searchsorted(self.erg_bins, erg, side='right') - 1

    @property
    def erg_bin_widths(self):
        return self.get_erg_bin_widths(adc=1)

    def get_erg_bin_widths(self, adc=1):
        assert adc > 0, "ADC starts at 1, not 0."
        return self._erg_bins[adc-1][1:] - self._erg_bins[adc-1][1:]

    def plot_spectrum(self, adc=1, ax=None, leg_name=None):
        if ax is None:
            plt.figure()
            ax = plt.gca()

        ax.errorbar(self.get_energies(adc), self.get_counts(adc), np.sqrt(self.get_counts(adc)), label=leg_name)
        if leg_name:
            ax.legend()


class MPANTList:
    def __init__(self, path, manual_mpa_path=None, dead_time_corr_window=20):
        """

        Args:
            path:
            manual_mpa_path:
            dead_time_corr_window: The width of convolution window for calculating % live time
        """
        path = Path(path)
        self.file_name = path.name
        self._energies = []
        self._times = []
        self._livetimes = []

        mpa_path = path.with_suffix(".mpa")

        if manual_mpa_path is None:
            try:
                self.mca = MPA(mpa_path)
            except FileNotFoundError:
                assert False, "Corresponding MPA file not found (used for erg calibration info, etc.). " \
                              "Use `manual_mpa_path` to provide MPA file with correct information. "
        else:
            self.mca = MPA(manual_mpa_path)

        prev_times = [0., 0., 0., 0.]
        tot_live_time = 0

        with open(path) as f:
            while f:
                line = f.readline()
                if not line:
                    break
                adc_index, ch, time = map(int, line.split())
                adc = adc_index + 1
                time = time*5E-8
                energy = self.mca.channel2erg(ch, adc)
                tot_live_time += time - prev_times[adc_index] - 1.25E-5
                prev_times[adc_index] = time
                try:
                    self._energies[adc_index].append(energy)
                    self._times[adc_index].append(time)
                    self._livetimes[adc_index].append(tot_live_time)
                except IndexError:
                    for i in range(len(self._times), adc_index + 1):
                        self._energies.append([])
                        self._times.append([])
                        self._livetimes.append([])
                    self._energies[adc_index].append(energy)
                    self._times[adc_index].append(time)
                    self._livetimes[adc_index].append(tot_live_time)

        self._percent_live: List[np.ndarray] = [np.ndarray([0]) for i in range(len(self._livetimes))]

        for adc_index in range(len(self._livetimes)):
            if not len(self._livetimes[adc_index]):
                continue

            kernel = np.ones(dead_time_corr_window//2)/(dead_time_corr_window//2)
            self._percent_live[adc_index] = \
                np.gradient(self._livetimes[adc_index]) / np.gradient(self._realtimes[adc_index])
            self._percent_live[adc_index] = np.convolve(self._percent_live[adc_index], kernel, mode='same')

            # correct edge effects
            self._percent_live[adc_index][0:dead_time_corr_window//2] = \
                np.median(self._percent_live[adc_index][dead_time_corr_window//2:dead_time_corr_window])
            # print(self._percent_live[adc_index][:20])

            self._percent_live[adc_index][-dead_time_corr_window//2:] = \
                np.median(self._percent_live[adc_index][-dead_time_corr_window: -dead_time_corr_window//2])

    @property
    def _realtimes(self):
        return self._times

    def erg_bin_index(self, erg, adc=1):
        assert adc > 0, "ADC starts at 1, not 0"
        return np.searchsorted(self.mca.get_erg_bins(adc=adc), erg, side='right') - 1

    def get_times(self, adc=1):
        assert adc > 0, "ADC starts at 1, not 0"
        return np.array(self._times[adc-1])

    @property
    def times(self):
        return self.get_times(adc=1)

    def get_energies(self, adc=1):
        assert adc > 0, "ADC starts at 1, not 0"
        return np.array(self._energies[adc-1])

    @property
    def energies(self):
        return self.get_energies(adc=1)

    def get_erg_bins(self, adc=1):
        return self.mca.get_erg_bins(adc=adc)

    @property
    def percent_live(self):
        return self.get_percent_live(adc=1)

    def get_percent_live(self, adc=1):
        assert adc > 0, "ADC starts at 1, not 0"
        return self._percent_live[adc-1]

    @property
    def livetimes(self):
        return self.get_livetimes(adc=1)

    def get_livetimes(self, adc=1):
        assert adc > 0, "ADC starts at 1, not 0"
        return self._livetimes[adc - 1]

    @property
    def erg_bins(self):
        return self.get_erg_bins(adc=1)

    @cached_property
    def erg_centers(self):
        return np.array([(b0 + b1) / 2 for b0, b1 in zip(self.erg_bins[:-1], self.erg_bins[1:])])

    @property
    def deadtime_corrs(self):
        return self.get_deadtime_corrs(adc=1)

    def get_deadtime_corrs(self, adc=1):
        assert adc > 0, "ADC starts at 1, not 0"
        a = self.get_percent_live(adc)
        return 1.0 / a

    def plot_count_rate(self, adc=None, bins=None):
        if bins is None:
            bins = 'auto'

        if adc is None:
            fig, axs = plt.subplots(2, 1, sharex='all')

            _times = np.concatenate(self._times)
        else:
            _times = np.concatenate(self._times[adc])
            axs = [plt.gca()]
        data_array = [self.get_times(adc)] if adc is not None else self._times

        for index, (times, ax) in enumerate(zip(data_array, axs)):
            y, bin_edges = np.histogram(times, bins)
            bin_widths = bin_edges[1:] - bin_edges[:-1]
            yerr = np.sqrt(bin_widths)

            y = y / bin_widths
            yerr /= bin_widths
            mpl_hist(bin_edges, y, np.sqrt(y), ax=ax, label=f'ADC {index+1}')
            ax.legend()
            ax.set_ylabel("count rate [Hz]")
            if index == 0:
                ax.set_title("MPANT count rates")
            else:
                ax.set_xlabel("Time [s]")

        plt.subplots_adjust(hspace=0)

    def plot_percent_live(self, adc=1, ax=None, **ax_kwargs):
        if ax is None:
            plt.figure()
            ax = plt.gca()
        ax.plot(self._realtimes[adc-1], self.get_percent_live(adc), **ax_kwargs)
        ax.set_xlabel("Real time [s]")
        ax.set_ylabel("Fraction of time det. is live")
        return ax


if __name__ == '__main__':
    l = MPANTList('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/list_files/10-100-1k-10k-100k.txt')
    # l = MPANTList('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/list_files/100000.txt')

    l.plot_percent_live(adc=2)
    l.plot_count_rate(bins=100)
    plt.show()