"""
ad-hoc efficiency calibration for IAC run.

steps:
    1. Run save_eff_points(data_dir_name)
        where,
            data_dir_name: path to folder containing SPE files (rel. to pwd/cal_data).
    2. Run save_points




"""
from __future__ import annotations
from JSB_tools import Nuclide
import numpy as np
from matplotlib import pyplot as plt
import re
from pathlib import Path
from JSB_tools.list_reader import MaestroListFile
from JSB_tools.spe_reader import SPEFile
from datetime import datetime
from cal_data.cal_sources import CalSource
from uncertainties import unumpy as unp
from typing import Dict
import pickle

top_level_data_path = Path.cwd()/'cal_data'  # where data is saved.


cal_source_serial_numbers = {'Na22': 129742, 'Mn54': 'J4-348',
               "Co57": "K4-895", "Cd109": '129757', 'Cs137':129792, 'Y88': 190607000}


#  ============================================
window_in_kev = 10  # Window width for counting peaks (No fitting is done here)
print_nuclide = 'Ag109_m1'  # Just for printing info about Nuclide.
plot_peaks = False  # For debugging.
#  ============================================

if __name__ == '__main__':
    print("Print Nuclide ===============================")
    _n = Nuclide.from_symbol(print_nuclide)
    print("\tDecay daughers: ", _n.decay_daughters)
    print(f"\t{_n.name} gammas: ")
    for g in _n.decay_gamma_lines:
        print('\t\t', g)
    date_of_cal = datetime(2021, 8, 17)
    print("End Print Nuclide ===============================")


def get_gammas(n: Nuclide, debug=False):
    """
    Return Nuclide.Gammaline instances for each nuclide. Edit this when adding more nuclides.
    Args:
        n:

    Returns:

    """
    n_parent = n
    # add special cases of daughter decays here.
    if n.name == 'Cd109':  # gamma comes from daughter nucleus
        n = Nuclide.from_symbol("Ag109_m1")
    elif n.name == 'Cs137':  # gamma comes from daughter nucleus
        n = Nuclide.from_symbol("Ba137_m1")
    else:
        n = n

    # add special cases of multiple nuclides here
    if n.name == 'Y88':
        out = n.decay_gamma_lines[:2]
    else:
        out = n.decay_gamma_lines[:1]

    if debug:
        print(f"Decay gamma lines for {n_parent}:")
        for g in n.decay_gamma_lines:
            print(f"\t{g}")
        print()
    return out


class EffCal:
    """
    #
        Attributes:
            points: entry for each source position. Entries are efficiencies. e.g. points[1][2] means 2nd position
                and third gamma line
            ergs: Energies of gamma lines under consideration.
    """
    instances: Dict[str, EffCal] = {}

    def effective_efficiencies(self):
        """
        Return weighted mean of efficiency over all source positions
        Returns:

        """
        points = [p for p in self.points if p]
        weights = np.ones(len(points))
        # weights[1:] += 1  # not sure about this...
        out = np.sum([eff_points*w for eff_points, w in zip(np.array(points), weights)], axis=0)
        out /= sum(weights)
        return out

    @classmethod
    def get_source(cls, name) -> EffCal:
        """Look up instance for a given source name."""
        return cls.instances[name]

    def __init__(self, nuclide_name, ergs):
        self.nuclide_name = nuclide_name

        self.points = [[], [], [], []]
        self.ergs = ergs
        EffCal.instances[nuclide_name] = self

    def add_point(self, true, meas, pos):
        """
        Add a gamma line from this source.
        Args:
            true:
            meas:
            pos:

        Returns:

        """
        self.points[pos].append(meas/true)


def save_eff_points(data_dir_name, debug=False):
    """

    Args:
        data_dir_name: Path to folder containing Spe files for the calibration. Format: e.g. Co-57-2, where 2 means the
            2nd position is being used. Position "0" corresponds to center of detector.
        debug:

    Returns:

    """
    data_path = top_level_data_path/data_dir_name
    bg_spe = SPEFile(data_path/'BG.Spe')
    back_ground_count_rates = bg_spe.counts / bg_spe.livetime
    plt.title("Background count rate")
    plt.plot(bg_spe.energies, unp.nominal_values(back_ground_count_rates))

    file_matches = []  # files that match SPE

    # sort files so that sources are processed together
    for f in sorted(data_path.iterdir()):
        if m := re.match("([A-Z][a-z]?-[0-9]+)-([0-9])\.Spe", f.name):
            file_matches.append((m, f))

    cashed_gammas = {}
    all_spe_center_files = {}  # all SPE files at center of detector (for debug only).

    for m, f in file_matches:
        nuclide = Nuclide.from_symbol(m.groups()[0])

        src_pos_i = int(m.groups()[1])  # Index corresponding to location of source for these data.
        cal_source = CalSource.get_source(cal_source_serial_numbers[nuclide.name])
        spe = SPEFile(f)

        if nuclide.name not in all_spe_center_files:
            all_spe_center_files[nuclide.name] = spe

        n_decays = cal_source.get_n_decays(spe.livetime, date_of_cal)

        spectrum_array = spe.counts - back_ground_count_rates*spe.livetime
        baseline = SPEFile.calc_background(spectrum_array)
        spectrum_array -= baseline

        try:
            cal_gammas = cashed_gammas[nuclide.name]
        except KeyError:
            cal_gammas = get_gammas(nuclide, debug)
            cashed_gammas[nuclide.name] = cal_gammas

        try:
            cal_dataclass = EffCal.get_source(nuclide.name)
        except KeyError:
            cal_dataclass = EffCal(nuclide.name, [g.erg.n for g in cal_gammas])

        for g in cal_gammas:
            n_gammas_true = n_decays * g.intensity
            ilow, ihigh = spe.erg_bin_index(g.erg.n-window_in_kev//2), spe.erg_bin_index(g.erg.n+window_in_kev//2)
            peak_counts = spectrum_array[ilow: ihigh]  # array of baseline subtracted counts in the neighborhood of peak
            meas_counts = sum(peak_counts)  # measured peak counts.
            cal_dataclass.add_point(n_gammas_true, meas_counts, src_pos_i)

            if src_pos_i == 0 and (plot_peaks or debug):
                plt.figure()
                plt.plot(spe.energies[ilow: ihigh], unp.nominal_values(peak_counts))
                plt.title(f"{nuclide.name}, center")

    plt.figure()

    eff_final_x = []
    eff_final_y = []
    eff_final_y_err = []
    for index_color, cal in enumerate(EffCal.instances.values()):  # loop through different sources.
        color = ['red', 'blue', 'green', 'black', 'purple', 'orange', 'yellow'][index_color]

        eff_final_x.extend(cal.ergs)
        eff_final_y.extend(unp.nominal_values(cal.effective_efficiencies()))
        eff_final_y_err.extend(unp.std_devs(cal.effective_efficiencies()))
        for index_marker, p in enumerate(cal.points):  # loop through different source positions
            p_err = [_.std_dev for _ in p]
            p = [_.n for _ in p]
            if p:  # plot all points for this nuclide at this position
                plt.errorbar(cal.ergs, p, yerr=p_err, label=f"{cal.nuclide_name}" if index_marker == 0 else None,
                             marker=[".", 'd', 'x', 'p'][index_marker],
                             color=color, ls='None')

    eff_final_x.insert(0, 5)
    eff_final_y.insert(0, 1E-5)
    eff_final_y_err.insert(0, eff_final_y_err[0]/3.)
    arg_sort = np.argsort(eff_final_x)
    eff_final_x = np.array(eff_final_x)[arg_sort]
    eff_final_y = np.array(eff_final_y)[arg_sort]
    eff_final_y_err = np.array(eff_final_y_err)[arg_sort]
    plt.legend()

    with open(data_path/'points.pickle', 'wb') as f:
        _d = {'x': eff_final_x, 'y': eff_final_y, 'y_err': eff_final_y_err}
        pickle.dump(_d, f)

    if debug:
        n = len(all_spe_center_files)
        rows = (n // 3 + (1 if n % 3 != 0 else 0))
        fig, axs = plt.subplots(rows, 3, figsize=(8,4), sharex='all')
        axs = axs.flatten()
        for index, (name, _spe) in enumerate(all_spe_center_files.items()):
            ax = axs[index]
            y = _spe.counts/_spe.livetime - back_ground_count_rates
            y_err = unp.std_devs(y)
            y = unp.nominal_values(y)
            ax.errorbar(_spe.energies, y, y_err, label=f"{name}")

            if index >= len(all_spe_center_files)-3:
                ax.set_xlabel("Energy [KeV]")
            ax.set_ylabel("Count rate [Hz]")
            ax.legend()


if __name__ == '__main__':
    save_eff_points("our_det/08_17", True)
    plt.show()
    # n = 5 // 2
    # rows = (n // 3 + (1 if n % 3 != 0 else 0))
    # fig, axs = plt.subplots(rows, 3)
    # axs = axs.flatten()
    # print(rows)
#

