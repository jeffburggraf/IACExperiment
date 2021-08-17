from pathlib import Path
import numpy as np
import re
from matplotlib import pyplot as plt
from datetime import datetime
from scipy.stats import norm


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
        self._counts = counts

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
        return self.energies(adc)[int(ch)]

    def live_time(self, adc=1):
        return self._live_times[adc-1]

    def real_time(self, adc=1):
        return self._real_times[adc-1]

    def counts(self, adc=1):
        return self._counts[adc-1]

    def energies(self, adc=1):
        return self._energies[adc-1]

    def get_erg_bins(self, adc=1):
        return self._erg_bins[adc-1]

    @property
    def erg_bins(self):
        return self.get_erg_bins(adc=1)

    def plot_spectrum(self, adc=1, ax=None, leg_name=None):
        if ax is None:
            plt.figure()
            ax = plt.gca()

        ax.errorbar(self.energies(adc), self.counts(adc), np.sqrt(self.counts(adc)), label=leg_name)
        if leg_name:
            ax.legend()


class MPANTList:
    # Todo: Maybe deadtime is simply tot_time-a*n_events? Test this.

    def __init__(self, path):
        path = Path(path)
        self.file_name = path.name
        self._energies = []
        self._times = []
        try:
            mpa_path = path.with_suffix(".mpa")
        except FileNotFoundError:
            assert False, f"MPA file corresponding to {path.name} not found!"

        self.mca = MPA(mpa_path)
        with open(path) as f:
            while f:
                line = f.readline()
                if not line:
                    break
                adc_index, ch, time = map(int, line.split())
                adc = adc_index + 1
                time = time*5E-8
                energy = self.mca.channel2erg(ch, adc)
                try:
                    self._energies[adc_index].append(energy)
                    self._times[adc_index].append(time)
                except IndexError:
                    try:
                        self._energies.insert(adc_index, [energy])
                        self._times.insert(adc_index, [time])
                    except IndexError:
                        self._energies.append([energy])
                        self._times.append([time])

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
    def erg_bins(self):
        return self.get_erg_bins(adc=1)

    def dead_time_corr(self, t, adc=1):
        return 1  # todo

    def plot_count_rate(self, adc=1):
        # todo
        pass


if __name__ == '__main__':
    # l = MPANTList('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/list_files/beamgun003.txt')
    # print(l.energies())
    bins = [1,2,3,5,6,7]
    print(np.searchsorted(bins, 3, side="right")-1)

