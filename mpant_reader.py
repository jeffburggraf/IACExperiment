import warnings
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
from typing import List, Dict
from uncertainties import umath


class MPA:
    def __init__(self, path, primary_adc_number=1, aux_adc_number=3):
        assert primary_adc_number > 0, '<=0 is not a valid ADC index'
        assert aux_adc_number > 0, '<=0 is not a valid ADC index'
        self.primary_adc_number = primary_adc_number
        self.aux_adc_number = aux_adc_number
        with open(path) as f:
            lines = f.readlines()
        header = lines[:lines.index("[DATA0,8192 ]\n")]
        adc_headers = {}
        for adc_number in range(1, 6):
            try:
                start_i = header.index(f"[ADC{adc_number}]\n")
                end_s = f"[ADC{adc_number+1}]\n"
                if end_s not in header:
                    adc_headers[adc_number] = header[start_i+1:]
                else:
                    end_i = header.index(end_s)
                    adc_headers[adc_number] = header[start_i+1: end_i]
            except ValueError:
                continue

        adc_headers_dict = {}
        for (adc_number, head) in adc_headers.items():
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
                adc_headers_dict[adc_number] = d

        counts = {adc_number: [] for adc_number in adc_headers_dict.keys()}

        adc_number = None
        for line in lines[len(header):]:
            if m := re.match("\[DATA([0-9]),", line):
                 adc_number = int(m.groups()[0]) + 1
            else:
                if adc_number is not None:
                    counts[adc_number].append(int(line))
                else:
                    warnings.warn("Expected a line like [DATA[0-9],[0-9]+ ] after header. Didnt see one. Beware.")

        self._counts = {adc_number: unp.uarray(c, np.sqrt(c)) for adc_number, c in counts.items()}

        channels = {adc_number: np.arange(0, adc_dict['range']) + 0.5 for adc_number, adc_dict in adc_headers_dict.items()}

        channels_bins = {adc_number: np.arange(0, adc_dict['range']+1) for adc_number, adc_dict in adc_headers_dict.items()}

        self._erg_coeffs = {adc_number: np.array([adc_dict['caloff'], adc_dict['calfact']])
                            for adc_number, adc_dict in adc_headers_dict.items()}

        self._energies = {adc_number: np.sum([c*chs**pow for pow, c in enumerate(self._erg_coeffs[adc_number])], axis=0)
                          for adc_number, chs in channels.items()}

        self._erg_bins = {adc_number: np.sum([c*ch_bins**pow for pow, c in enumerate(self._erg_coeffs[adc_number])],
                          axis=0) for adc_number, ch_bins in channels_bins.items()}

        self._live_times = {adc_number: adc_dict['livetime'] for adc_number, adc_dict in adc_headers_dict.items()}
        self._real_times = {adc_number: adc_dict['realtime'] for adc_number, adc_dict in adc_headers_dict.items()}
        # self._real_times = [adc_dict['runtime'] for adc_dict in adc_headers_dict.values()]

        for line in header:
            if m := re.match("REPORT-FILE from (.+) written", line):
                s = m.groups()[0].strip()
                self.system_start_time = datetime.strptime(s, "%m/%d/%Y %H:%M:%S")
                break

    @cached_property
    def valid_adcs(self):
        return self._live_times.keys()

    def __get_adc_index__(self, adc):
        if adc is None:
            return self.primary_adc_number
        else:
            assert adc > 0,  "ADC starts at 1, not 0."
            return adc

    def channel2erg(self, ch, adc=1):
        # todo: Is this right? is there a "zero" channel? Compare to MPANT.
        return self.get_energies(adc)[int(ch)]

    @property
    def livetime(self):
        return self.get_livetime()

    def get_livetime(self, adc=None):
        return self._live_times[self.__get_adc_index__(adc)]

    @property
    def realtime(self):
        return self.get_realtime()

    def get_realtime(self, adc=None):
        if adc is None:
            adc = self.primary_adc_number
        else:
            assert adc > 0, "ADC starts at 1, not 0."
        return self._real_times[adc]

    @property
    def counts(self, adc=None, nominal=False):
        return self.get_counts(self.__get_adc_index__(adc), nominal=nominal)

    def get_counts(self, adc=None, nominal=False):
        adc = self.__get_adc_index__(adc)
        if nominal:
            return unp.nominal_values(self._counts[adc])
        else:
            return self._counts[adc]

    @property
    def energies(self, adc=None):
        return self.get_energies(self.__get_adc_index__(adc))

    def get_energies(self, adc=None):
        return self._energies[self.__get_adc_index__(adc)]

    @property
    def erg_bins(self):
        return self.get_erg_bins(adc=self.primary_adc_number)

    def get_erg_bins(self, adc=None):
        return self._erg_bins[self.__get_adc_index__(adc)]

    def erg_bin_index(self, erg):
        if hasattr(erg, '__iter__'):
            return np.array([self.erg_bin_index(e) for e in erg])
        if isinstance(erg, AffineScalarFunc):
            erg = erg.n
        return np.searchsorted(self.erg_bins, erg, side='right') - 1

    @property
    def erg_bin_widths(self):
        return self.get_erg_bin_widths(adc=self.primary_adc_number)

    def get_erg_bin_widths(self, adc=None):
        adc = self.__get_adc_index__(adc)
        return self._erg_bins[adc][1:] - self._erg_bins[adc][1:]

    def get_energies_in_range(self, erg_min, erg_max, adc=None):
        """
        Integrate events over energy.
        Args:
            erg_min:
            erg_max:
            adc:

        Returns:

        """
        adc = self.__get_adc_index__(adc)
        return self.get_energies(adc)[np.where((self.get_energies(adc) >= erg_min) & (self.get_energies(adc) <= erg_max))]

    def plot_spectrum(self, adc=None, ax=None, leg_name=None, erg_min=None, erg_max=None,
                      **mpl_kwargs):
        adc = self.__get_adc_index__(adc)
        if ax is None:
            plt.figure()
            ax = plt.gca()
        if erg_min is not None:
            min_index = self.erg_bin_index(erg_min)
        else:
            min_index = 0
        if erg_max is not None:
            max_index = self.erg_bin_index(erg_max)
        else:
            max_index = len(self.erg_bins)-1
        # ax.errorbar(self.get_energies(adc)[min_index: max_index],
        #             self.get_counts(adc, nominal=True)[min_index: max_index],
        #             unp.std_devs(self.get_counts(adc))[min_index: max_index],
        #             label=leg_name)
        mpl_hist(self.erg_bins[min_index: max_index+1], self.get_counts(adc)[min_index: max_index],
                 ax=ax, label=leg_name, **mpl_kwargs)
        ax.set_xlabel('energy [KeV]')
        ax.set_ylabel('counts')
        if leg_name:
            ax.legend()


class MPANTList:
    def __init__(self, path, manual_mpa_path=None, primary_adc_number=1, aux_adc_number=3, dead_time_corr_window=20):
        """

        Args:
            path:
            manual_mpa_path:
            dead_time_corr_window: The width of convolution window for calculating % live time
        """
        path = Path(path)
        self.file_name = path.name
        self._energies = {}
        self._times = {}
        self._livetimes = {}

        mpa_path = path.with_suffix(".mpa")
        self.primary_adc_number = primary_adc_number
        self.aux_adc_number = aux_adc_number
        if manual_mpa_path is None:
            try:
                self.mca = MPA(mpa_path, primary_adc_number=primary_adc_number, aux_adc_number=aux_adc_number)
            except FileNotFoundError:
                raise FileNotFoundError("Corresponding MPA file not found (used for erg calibration info, etc.). "
                                        "Use `manual_mpa_path` to provide MPA file with correct information. ")
        else:
            self.mca = MPA(manual_mpa_path)

        prev_times = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
        tot_live_time: float = 0

        with open(path) as f:
            while f:
                line = f.readline()
                if not line:
                    break
                _adc_index_, ch, time = map(int, line.split())
                adc_number = _adc_index_ + 1
                time = time*5E-8
                energy = self.mca.channel2erg(ch, adc_number)
                tot_live_time += time - prev_times[adc_number] - 1.25E-5
                prev_times[adc_number] = time
                try:
                    self._energies[adc_number].append(energy)
                    self._times[adc_number].append(time)
                    self._livetimes[adc_number].append(tot_live_time)
                except KeyError:
                    self._energies[adc_number] = [energy]
                    self._times[adc_number] = [time]
                    self._livetimes[adc_number] = [tot_live_time]

        self._percent_live: Dict[int, np.ndarray] = {adc_number: np.ndarray([0])
                                                     for adc_number in self._livetimes.keys()}

        for adc_number in self._livetimes.keys():
            if not len(self._livetimes[adc_number]):
                continue

            kernel = np.ones(dead_time_corr_window//2)/(dead_time_corr_window//2)
            self._percent_live[adc_number] = \
                np.gradient(self._livetimes[adc_number]) / np.gradient(self._realtimes[adc_number])
            self._percent_live[adc_number] = np.convolve(self._percent_live[adc_number], kernel, mode='same')

            # correct edge effects
            self._percent_live[adc_number][0:dead_time_corr_window//2] = \
                np.median(self._percent_live[adc_number][dead_time_corr_window//2:dead_time_corr_window])
            # print(self._percent_live[adc_index][:20])

            self._percent_live[adc_number][-dead_time_corr_window//2:] = \
                np.median(self._percent_live[adc_number][-dead_time_corr_window: -dead_time_corr_window//2])

    @cached_property
    def valid_adcs(self):
        return self._times.keys()

    def __get_adc_index__(self, adc) -> int:
        if adc is None:
            adc = self.primary_adc_number
        else:
            assert adc > 0, "ADC starts at 1, not 0"
        assert adc in self.valid_adcs, f"Invalid ADC number, {adc}. The following are available:\n{self.valid_adcs}"
        return adc

    @property
    def _realtimes(self):
        return self._times

    def erg_bin_index(self, erg, adc=None):
        adc = self.__get_adc_index__(adc)
        return np.searchsorted(self.mca.get_erg_bins(adc=adc), erg, side='right') - 1

    def get_times(self, adc=None):
        adc = self.__get_adc_index__(adc)
        return np.array(self._times[adc])

    @property
    def times(self):
        return self.get_times(adc=None)

    def get_energies(self, adc=None):
        adc = self.__get_adc_index__(adc)
        return np.array(self._energies[adc])

    @property
    def energies(self):
        return self.get_energies(adc=None)

    def get_erg_bins(self, adc=None):
        return self.mca.get_erg_bins(adc=adc)

    @property
    def percent_live(self):
        return self.get_percent_live(adc=None)

    def get_percent_live(self, adc=None):
        adc = self.__get_adc_index__(adc)
        return self._percent_live[adc]

    @property
    def livetimes(self):
        return self.get_livetimes(adc=None)

    def get_livetimes(self, adc=None):
        adc = self.__get_adc_index__(adc)
        return self._livetimes[adc]

    @property
    def erg_bins(self):
        return self.get_erg_bins(adc=None)

    @cached_property
    def erg_centers(self):
        return np.array([(b0 + b1) / 2 for b0, b1 in zip(self.erg_bins[:-1], self.erg_bins[1:])])

    @property
    def deadtime_corrs(self):
        return self.get_deadtime_corrs(adc=None)

    def get_deadtime_corrs(self, adc=None):
        adc = self.__get_adc_index__(adc)
        a = self.get_percent_live(adc)
        return 1.0 / a

    def plot_count_rate(self, adc="All", bins=None):
        if bins is None:
            bins = 'auto'

        if adc == 'All':
            fig, axs = plt.subplots(2, 1, sharex='all')
            _times = np.concatenate(tuple(self._times.values()))
            adcs, data_array = self._times.keys(), self._times.values()
        else:
            adc = self.__get_adc_index__(adc)
            _times = np.concatenate(self._times[adc])
            axs = [plt.gca()]
            data_array = [self.get_times(adc)]
            adcs = [adc]
        # data_array = [self.get_times(adc)] if adc is not None else tuple(self._times.values())

        for index, (adc, times, ax)in enumerate(zip(adcs, data_array, axs)):
            y, bin_edges = np.histogram(times, bins)
            bin_widths = bin_edges[1:] - bin_edges[:-1]
            yerr = np.sqrt(bin_widths)

            y = y / bin_widths
            yerr /= bin_widths
            mpl_hist(bin_edges, y, np.sqrt(y), ax=ax, label=f'ADC {adc}')
            ax.legend()
            ax.set_ylabel("count rate [Hz]")
            if index == 0:
                ax.set_title("MPANT count rates")
            else:
                ax.set_xlabel("Time [s]")

        plt.subplots_adjust(hspace=0)

    def plot_percent_live(self, adc=None, ax=None, **ax_kwargs):
        adc = self.__get_adc_index__(adc)
        if ax is None:
            plt.figure()
            ax = plt.gca()
        ax.plot(self._realtimes[adc], self.get_percent_live(adc), **ax_kwargs)
        ax.set_xlabel("Real time [s]")
        ax.set_ylabel("Fraction of time det. is live")
        return ax


if __name__ == '__main__':
    l = MPANTList('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/list_files/jeff002.txt', )
    print(l.valid_adcs)
    # l.plot_percent_live(adc=3)
    plt.plot(l.get_livetimes(3), l.get_times(3))
    plt.show()
    # mpa = MPA('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/list_files/beamgun003.mpa')
    # # mpa = MPA('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/list_files/jeff002.mpa')
    # mpa.plot_spectrum(adc=1)
    # plt.show()

    # l.plot_percent_live(adc=2)
    # l.plot_count_rate(bins=100)
    # plt.show()