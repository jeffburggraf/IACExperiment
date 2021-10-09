from pathlib import Path
import numpy as np
from JSB_tools.spe_reader import SPEFile
from mpant_reader import MPA
from matplotlib import pyplot as plt
import uncertainties.unumpy as unp
from typing import Union
from JSB_tools import Nuclide, mpl_hist, shade_plot
from cal_data.cal_sources import CalSource
data_dir = Path(__file__).parent.parent/'exp_data'
from JSB_tools.regression import LogPolyFit

sources = {"Eu152": CalSource.get_source('129753'), "Co57": CalSource.get_source('k4-896'),
           'Cs137': CalSource.get_source(129792), 'Y88': CalSource.get_source(190607000)}


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
            self.spec: SPEFile = SPEFile(path)
            if cal_spe is not None:
                self.spec.set_energy_cal(*cal_spe.erg_calibration)
        else:
            self.spec: SPEFile = MPA(path)
        self.counts = self.spec.counts
        self.counts = self.spec.get_baseline_removed()
        self.rates = self.counts/self.spec.livetime
        self.energies = self.spec.energies
        self.true_decays = self.source.get_n_decays(self.spec.realtime, self.spec.system_start_time)

        _, label_eff = self.calc_efficiency(0)
        self.label = f'{nuclide};{"LLNL" if self.llnl_det else "IAC"};{self.day}-' \
                     f'{self.spec.system_start_time.time()};pos={self.position};eff={label_eff:.2e}'

    def plot_rate(self, ax=None, alpha=0.7, density=True):
        if ax is None:
            plt.figure()
            ax = plt.gca()
        if density:
            y = self.rates / self.spec.erg_bin_widths
        else:
            y = self.rates
        mpl_hist(self.spec.erg_bins, y, ax=ax, alpha=alpha, label=self.label)
        ax.set_xlabel('energy [KeV]')
        if not density:
            ax.set_ylabel('count rate [Hz]')
        else:
            ax.set_ylabel('Count rate [Hz/KeV]')
        ax.legend()

        return ax

    def calc_efficiency(self, gamma_index_or_energy: Union[float, int] = 0, window_kev=None, debug_plot=False):
        if window_kev is None:
            window_kev = 7 if self.llnl_det else 20
        if isinstance(gamma_index_or_energy, int):
            gamma_index = gamma_index_or_energy
            if self.nuclide.name != 'Cs137':
                gamma_line = self.nuclide.decay_gamma_lines[gamma_index]
            else:
                gamma_line = Nuclide.from_symbol('Ba137_m1').decay_gamma_lines[gamma_index]
        else:
            gamma_line = self.nuclide.get_gamma_nearest(gamma_index_or_energy)

        erg = gamma_line.erg.n
        # min_index = self.spec.__erg_index__(erg-window_kev/2)
        # max_index = self.spec.__erg_index__(erg+window_kev/2)
        # tot_events = np.sum(self.counts[min_index: max_index+1])
        tot_events = np.sum(self.spec.get_counts(erg-window_kev/2, erg+window_kev/2, remove_baseline=True,
                                                 deadtime_corr=True))
        true_events = self.true_decays*gamma_line.intensity
        out = tot_events/true_events
        if debug_plot is not None:
            ax = self.spec.plot_erg_spectrum(erg_min=erg-10*window_kev, erg_max=erg+10*window_kev, remove_baseline=True,
                                             make_rate=True, make_density=True)
            ax.set_title(f"{self.nuclide.name}: {gamma_line.erg} KeV\neff:{out}")
            shade_plot(ax, [erg-window_kev/2, erg+window_kev/2], label=f'window; n_events: {tot_events}')
            ax.legend()
        return erg, out

    def print_gamma_lines(self):
        for g in self.nuclide.decay_gamma_lines:
            print(g)


if __name__ == '__main__':
    """
    Notes:
    LLNL efficiency is pretty consistent (+/- 5%). Slight decrease in erg resolution on Friday, but that's it.
    
    """
    llnl_det = False
    X = []
    Y = []
    from JSB_tools.nuke_data_tools.gamma_coince import Levels
    Levels('Eu152').print_coinc(probability_cut_off=0.01)
    levels = Levels('Eu152')
    if llnl_det:
        eu_ergs = [295.939, 344.279, 964.059, 1112.081]

        c1 = CalData('wednesday', True, evening=False, nuclide='Co57')
        c2 = CalData('wednesday', True, evening=False, nuclide='Cs137')
        c2.plot_rate()
        c3 = CalData('wednesday', True, evening=True, nuclide='Eu152')
    else:
        c1 = CalData('thursday', False, evening=False, nuclide='Co57')
        c2 = CalData('thursday', False, evening=False, nuclide='Cs137')
        c2.plot_rate()
        c3 = CalData('thursday', False, evening=False, nuclide='Eu152')
        eu_ergs = [121.782, 344.279, 964.059, 1085., 1112.081,  1408.020]

    fig, ax = plt.subplots()
    for c in [c1, c3]:
        xs, ys = [], []
        if c.nuclide.name == 'Eu152':
            window = None if llnl_det else 40
            for e in eu_ergs:
                x, y, = c.calc_efficiency(gamma_index_or_energy=e, debug_plot=True, window_kev=25)
                xs.append(x)
                ys.append(y)
        else:
            for i in range(2):
                try:
                    x, y = c.calc_efficiency(gamma_index_or_energy=i, debug_plot=True)
                    if y.std_dev/y.n<0.25:
                        xs.append(x)
                        ys.append(y)
                except IndexError:
                    break
        X.extend(xs)
        Y.extend(ys)
        ax.errorbar(xs, unp.nominal_values(ys), unp.std_devs(ys), ls='None', marker='o', label=c.nuclide.name)

    srt = np.argsort(X)
    _X = np.array(X)[srt]
    _Y = np.array(Y)[srt]
    for e1 in eu_ergs:
        for e2, p in levels.get_coince(e1):
            e2_eff = np.interp(e2, _X, unp.nominal_values(_Y))
            s_map = np.isclose(e2 + e1, _X, atol=0.5)
            if any(s_map):
                print("Sum corr for ", _X[np.argmax(s_map)], p * e2_eff)
            if e2_eff*p>0.005:
                print('sub corr for ', e1, e2_eff*p)
    ax.legend()
    ax = c3.plot_rate()
    X.append(20)
    Y.append(Y[0]*1E-10)
    fit = LogPolyFit(X, Y, fix_coeffs=[0, 1, 2], order=3)
    fit.plot_fit()

    plt.show()