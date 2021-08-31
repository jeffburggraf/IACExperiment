from pathlib import Path
import numpy as np
from JSB_tools.spe_reader import SPEFile
from mpant_reader import MPA
from matplotlib import pyplot as plt
import uncertainties.unumpy as unp
from typing import Union
from JSB_tools import Nuclide
from cal_data.cal_sources import CalSource
data_dir = Path(__file__).parent.parent/'exp_data'

sources = {"Eu152": CalSource.get_source('129753')}


class CalData:
    def __init__(self, day: str, llnl_det=True, evening=False, nuclide='Eu152', position=0):
        self.day = day.lower()
        self.source: CalSource = sources[nuclide]
        self.llnl_det = llnl_det
        self.evening = evening
        self.nuclide = Nuclide.from_symbol(nuclide)
        self.position = position
        assert self.day in ['tuesday', 'wednesday', 'thursday', 'friday']
        cal_spe = None
        if llnl_det:
            dir_name = Path(self.day)/f'EffCal{self.day[0].upper()}{self.day[1:]}'
            suffix = '.Spe'
            if day == 'tuesday' and not evening:
                cal_spe = SPEFile(data_dir/dir_name/'use_for_cal.spe')
        else:
            dir_name = Path(self.day)/'MCA'/f'EffCal{self.day[0].upper()}{self.day[1:]}'
            suffix = '.mpa'

        if evening:
            dir_name /= "EndOfDay"
        dir_name /= f"{nuclide}-{position}.Spe"
        path = data_dir/dir_name
        path = path.with_suffix(suffix)
        if not path.exists():
            raise FileNotFoundError(path)
        self.spec: Union[SPEFile, MPA]
        if llnl_det:
            self.spec = SPEFile(path)
            if cal_spe is not None:
                self.spec.set_erg_cal(*cal_spe.erg_fit)
        else:
            self.spec = MPA(path)
        self.counts = self.spec.counts
        self.counts -= SPEFile.calc_background(self.counts)
        self.rates = self.counts/self.spec.livetime
        self.energies = self.spec.energies

        self.label = f'{nuclide};{"LLNL" if self.llnl_det else "IAC"};{self.day}-' \
                     f'{self.spec.system_start_time.time()};pos={self.position}'
    
    def plot_rate(self, ax=None, alpha=0.7):
        if ax is None:
            plt.figure()
            ax = plt.gca()
        ax.errorbar(self.energies, unp.nominal_values(self.rates), unp.std_devs(self.rates), label=self.label,
                    ds='steps-post', alpha=alpha)
        ax.legend()
        ax.set_xlabel('energy [KeV]')
        ax.set_ylabel('count rate [Hz]')
        return ax

    def calc_efficiency(self, gamma_index=0, window_kev=5):
        gamma_line = self.nuclide.decay_gamma_lines[gamma_index]
        erg = self.nuclide.decay_gamma_lines[gamma_index].erg.n
        min_index = self.spec.erg_bin_index(erg-window_kev//2)
        max_index = self.spec.erg_bin_index(erg+window_kev//2)
        tot_events = np.sum(self.counts[min_index: max_index+1])
        true_events = self.source.get_n_decays(self.spec.realtime, self.spec.system_start_time)*gamma_line.intensity
        return erg, tot_events/true_events

    def print_gamma_lines(self):
        for g in self.nuclide.decay_gamma_lines:
            print(g)


if __name__ == '__main__':
    c1 = CalData('friday', evening=False, llnl_det=True, nuclide='Eu152')
    c2 = CalData('tuesday', evening=False, llnl_det=True, nuclide='Eu152')
    c3 = CalData('wednesday', evening=True, llnl_det=True, nuclide='Eu152', position=0)
    c4 = CalData('wednesday', evening=False, llnl_det=False, nuclide='Eu152', position=0)
    c1.print_gamma_lines()

    ax = c1.plot_rate()
    c2.plot_rate(ax)
    c3.plot_rate(ax)
    c4.plot_rate(ax)
    for i in range(0, 6):
        print(c2.calc_efficiency(i))
    plt.show()