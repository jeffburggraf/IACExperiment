import warnings
from pathlib import Path
import numpy as np
import re
from matplotlib import pyplot as plt
from datetime import datetime
from functools import cached_property
from JSB_tools import mpl_hist
from uncertainties import unumpy as unp
from typing import Dict
from JSB_tools.spe_reader import SPEFile


class MPA(SPEFile):
    """
    Builds an SPEFile class from the MPA data and returns an instance of that class.
    """
    @classmethod
    def from_pickle(cls, f_path) -> SPEFile:
        path = (f_path.parent / f'{f_path.with_suffix("").name}_mpa').with_suffix('.Spe')
        return SPEFile.from_pickle(path)

    def __new__(cls, *args, **kwargs) -> SPEFile:
        self = super(MPA, cls).__new__(cls)
        self.__init__(*args, **kwargs)
        # path = (self.path.parent / f'{self.path.with_suffix("").name}_mpa').with_suffix('.Spe')
        spe = cls.build(path=self.path,
                        counts=self.mpa_counts[self.primary_adc_number],
                        erg_calibration=self.mpa_erg_coeffs[self.primary_adc_number],
                        livetime=self.mpa_live_times[self.primary_adc_number],
                        realtime=self.mpa_real_times[self.primary_adc_number],
                        channels=self.mpa_channels[self.primary_adc_number],
                        system_start_time=self.system_start_time, eff_path=self.eff_path
                        )
        return spe

    def __init__(self, mca_path, eff_path=None, primary_adc_number=1, aux_adc_number=3):
        assert primary_adc_number > 0, '<=0 is not a valid ADC index'
        assert aux_adc_number > 0, '<=0 is not a valid ADC index'
        mca_path = Path(mca_path)
        self.eff_path = Path(eff_path) if eff_path is not None else None
        self.path = (mca_path.parent / f'{mca_path.with_suffix("").name}_mpa').with_suffix('.Spe')
        self.primary_adc_number = primary_adc_number
        self.aux_adc_number = aux_adc_number
        with open(mca_path) as f:
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
                    warnings.warn("Expected a line that would match the following regex after the header: "
                                  "'[DATA[0-9],[0-9]+ ]'. Didn't see one, so beware.")

        self.mpa_counts = {adc_number: unp.uarray(c, np.sqrt(c)) for adc_number, c in counts.items()}

        self.mpa_channels = {adc_number: np.arange(0, adc_dict['range'], dtype=int) for adc_number, adc_dict in adc_headers_dict.items()}

        channels_bins = {adc_number: np.arange(0, adc_dict['range']+1)-0.5 for adc_number, adc_dict in adc_headers_dict.items()}

        self.mpa_erg_coeffs = {adc_number: np.array([adc_dict['caloff'], adc_dict['calfact']])
                            for adc_number, adc_dict in adc_headers_dict.items()}

        self.mpa_energies = {adc_number: np.sum([c*chs**pow for pow, c in enumerate(self.mpa_erg_coeffs[adc_number])], axis=0)
                          for adc_number, chs in self.mpa_channels.items()}

        self.mpa_erg_bins = {adc_number: np.sum([c*ch_bins**pow for pow, c in enumerate(self.mpa_erg_coeffs[adc_number])],
                          axis=0) for adc_number, ch_bins in channels_bins.items()}

        self.mpa_live_times = {adc_number: adc_dict['livetime'] for adc_number, adc_dict in adc_headers_dict.items()}
        self.mpa_real_times = {adc_number: adc_dict['realtime'] for adc_number, adc_dict in adc_headers_dict.items()}

        for line in header:
            if m := re.match("REPORT-FILE from (.+) written", line):
                s = m.groups()[0].strip()
                self.system_start_time = datetime.strptime(s, "%m/%d/%Y %H:%M:%S")
                break


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
    pass