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
from JSB_tools.regression import LogPolyFit

sources = {"Eu152": CalSource.get_source('129753'), "Co57": CalSource.get_source('k4-896')}


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
                self.spec.set_erg_cal(*cal_spe.erg_calibration)
        else:
            self.spec = MPA(path)
        self.counts = self.spec.counts
        self.counts -= SPEFile.calc_background(self.counts)
        self.rates = self.counts/self.spec.livetime
        self.energies = self.spec.energies

        self.label = f'{nuclide};{"LLNL" if self.llnl_det else "IAC"};{self.day}-' \
                     f'{self.spec.system_start_time.time()};pos={self.position}'
        self.true_decays = self.source.get_n_decays(self.spec.realtime, self.spec.system_start_time)

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

    def calc_efficiency(self, gamma_index=0, window_kev=None):
        if window_kev is None:
            window_kev = 5 if self.llnl_det else 15
        gamma_line = self.nuclide.decay_gamma_lines[gamma_index]
        erg = self.nuclide.decay_gamma_lines[gamma_index].erg.n
        min_index = self.spec.erg_bin_index(erg-window_kev//2)
        max_index = self.spec.erg_bin_index(erg+window_kev//2)
        tot_events = np.sum(self.counts[min_index: max_index+1])
        true_events = self.true_decays*gamma_line.intensity
        return erg, tot_events/true_events

    def print_gamma_lines(self):
        for g in self.nuclide.decay_gamma_lines:
            print(g)


if __name__ == '__main__':
    """
    Notes:
    LLNL Thurday: Same throughout the day by within a few %
    
    """
    eullnl = CalData('friday', evening=False, llnl_det=True, nuclide='Eu152')
    eullnl2 = CalData('thursday', evening=False, llnl_det=True, nuclide='Eu152')
    eullnl3 = CalData('thursday', evening=True, llnl_det=True, nuclide='Eu152')
    eullnl4 = CalData('wednesday', evening=True, llnl_det=True, nuclide='Eu152')
    eullnl5 = CalData('tuesday', evening=False, llnl_det=True, nuclide='Eu152')
    collnl = CalData('friday', evening=False, llnl_det=True, nuclide='Co57')
    eu_iac1 = CalData('wednesday', evening=False, llnl_det=False, nuclide='Eu152', position=0)
    eu_iac2 = CalData('friday', evening=False, llnl_det=False, nuclide='Eu152', position=0)
    eu_iac3 = CalData('thursday', evening=False, llnl_det=False, nuclide='Eu152', position=0)

    ax = eullnl.plot_rate()
    # collnl.plot_rate(ax)
    # eu_iac_w_morn.plot_rate(ax)
    # eu_iac_th_morn.plot_rate(ax)
    eullnl2.plot_rate(ax)
    eullnl3.plot_rate(ax)
    eullnl4.plot_rate(ax)
    eullnl5.plot_rate(ax)
    cal_points_x = [30]
    cal_points_y = [1E-6]

    for t in [eullnl, eullnl2, eullnl3, eullnl4, eullnl5, eu_iac1, eu_iac2, eu_iac3]:
        print(t.label, t.calc_efficiency())

    for i in range(1, 6):
        erg, eff = eullnl.calc_efficiency(i)
        cal_points_x.append(erg)
        cal_points_y.append(eff)
    for i in range(2):
        erg, eff = collnl.calc_efficiency(i)
        cal_points_x.append(erg)
        cal_points_y.append(eff)
    plt.figure()
    plt.errorbar(cal_points_x, unp.nominal_values(cal_points_y), unp.std_devs(cal_points_y), ls='None', marker='o')
    l1 = LogPolyFit(cal_points_x, cal_points_y, order=4, fix_coeffs=[0, 3])
    print(l1.fit_result.fit_report())
    l1.plot_fit()


    plt.show()