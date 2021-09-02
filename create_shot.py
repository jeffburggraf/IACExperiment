from __future__ import annotations
import pickle
from pathlib import Path
import ROOT
import pylab as p

from JSB_tools import ROOT_loop, TBrowser, mpl_hist
import numpy as np
from JSB_tools.list_reader import MaestroListFile
from matplotlib import pyplot as plt
from typing import Union, List
from mpant_reader import MPANTList
import dill
import marshal
from eff_calibration_setup import top_level_data_path
from lmfit.model import ModelResult
from JSB_tools import ProgressReport
from functools import cached_property
from uncertainties import unumpy as unp
from uncertainties import ufloat
list_file_dir = Path.cwd()/'exp_data'/'list_files'

# import matplotlib
# matplotlib.use('TkAgg')




class IACShot:
    """
    Nicely package up a IAC shot.

    To re-tool this for PHELIX, just change the get_beam_times function, and allow for both list file arguments
        to be maestro.

    Attributes:
        _list_of_data: 2D list of energies and times. Shape: (# erg channels, number of events in channel [varies]))
    """
    root_files = []
    save_path = Path.cwd()/'exp_data'/'packaged_shots'

    @classmethod
    def gen_fake_data(cls):
        from JSB_tools import Nuclide
        xe = Nuclide.from_symbol('Xe139')

        mean_transport_time = 10
        sigma_transport_time = 4
        sig_2_noise = 10
        n_samples = 200000

        gammas = xe.decay_gamma_lines
        gamma_ergs = [g.erg.n for g in gammas]
        probs = np.array([g.intensity.n for g in gammas])
        probs /= sum(probs)

        def gen_decay_ergs():
            return np.random.choice(gamma_ergs, size=n_samples, p=probs)

        def time_reject(t):
            sampled_t = np.random.normal(t, sigma_transport_time, 1)
            return sampled_t > mean_transport_time

        out = cls.__new__(cls)

        actual_ergs = gen_decay_ergs()
        meas_ergs = []
        actual_times = np.random.exponential(1.0/xe.decay_rate.n, size=n_samples)
        max_erg = max(actual_ergs)
        bg_times = np.random.uniform(min(actual_times), max(actual_times), size=n_samples)

        meas_times = []
        noise_prob = 1.0/(1+sig_2_noise)
        out.erg_bins = np.arange(0, max_erg+0.25, 0.25)

        out._list_of_data = []
        for i in range(len(out.erg_bins)-1):
            out._list_of_data.append([])

        for index, (t_meas, t_bg, erg) in enumerate(zip(actual_times, bg_times, actual_ergs)):
            if time_reject(t_meas):
                erg = np.random.normal(erg, 3)
                if erg>=max_erg:
                    continue
                meas_ergs.append(erg)
                meas_times.append(t_meas)
                out._list_of_data[np.searchsorted(out.erg_bins, erg, side='right') - 1].append(t_meas)

                if np.random.random() < noise_prob:
                    erg = np.random.uniform(min(gamma_ergs), max(gamma_ergs))
                    if erg >= max_erg:
                        continue
                    meas_ergs.append(erg)
                    meas_times.append(t_bg)
                    out._list_of_data[np.searchsorted(out.erg_bins, erg, side='right')-1].append(t_bg)
            if index%10000 == 0:
                print(index, erg)
        out.eff_fit = lambda x: 1
        out.event_energies = meas_ergs
        out.event_times = meas_times
        out.beam_center_time = 0
        out.beam_abs_times = [0]
        out.beam_duration = 0.1
        root_path, pickle_path, marshal_path = out.__get_paths__('fake')
        out.pickle(pickle_path, marshal_path)
        root_file = ROOT.TFile(str(root_path), 'recreate')
        IACShot.root_files.append(root_file)
        out._build_tree()

        return out

    @staticmethod
    def get_beam_times(list_file_maestro: MaestroListFile, list_file_mpant: MPANTList):
        beam_abs_times_maestro = list_file_maestro.realtimes[np.where(list_file_maestro.sample_ready_state == 1)]
        beam_abs_times_mpant = list_file_mpant.get_times(adc=list_file_mpant.aux_adc_number)

        beam_on_center_time_maestro = np.median(beam_abs_times_maestro)
        beam_on_center_time_mpant = np.median(beam_abs_times_mpant)
        beam_duration_maestro = 2 * np.percentile(beam_abs_times_maestro - beam_on_center_time_maestro, 98)
        beam_duration_mpant = 2 * np.percentile(beam_abs_times_mpant - beam_on_center_time_mpant, 98)
        return (beam_on_center_time_maestro, beam_on_center_time_mpant), \
               (beam_duration_maestro, beam_duration_mpant), \
               (beam_abs_times_maestro, beam_abs_times_mpant)

    @staticmethod
    def __get_paths__(name):
        root_path = (IACShot.save_path / name).with_suffix(".root")
        pickle_path = (IACShot.save_path / name).with_suffix(".pickle")
        marshal_path = (IACShot.save_path / name).with_suffix(".marshal")
        return root_path, pickle_path, marshal_path

    def _build_tree(self):
        br_time = np.array([0.], dtype=np.float32)
        br_energy = np.array([0.], dtype=np.float32)
        br_eff = np.array([0.], dtype=np.float32)
        br_dead_corr = np.array([0.], dtype=np.float32)

        self.trees = ROOT.TTree('tree1', 'tree1'), ROOT.TTree('tree2', 'tree2')

        self.trees[0].Branch("t", br_time, 't/F')
        self.trees[0].Branch("erg", br_energy, 'erg/F')
        self.trees[0].Branch("eff", br_eff, 'eff/F')
        self.trees[0].Branch('dead_corr', br_dead_corr, 'dead_corr/F')

        self.trees[1].Branch("t", br_time, 't/F')
        self.trees[1].Branch("erg", br_energy, 'erg/F')
        self.trees[1].Branch("eff", br_eff, 'eff/F')
        self.trees[1].Branch('dead_corr', br_dead_corr, 'dead_corr/F')

        eff_lists = [self.eff_fits[i].eval(x=self.event_energies[i]) for i in range(2)]

        prog = ProgressReport(len(self.event_energies[0])+len(self.event_energies[1]))
        prog_index = 0
        
        for det_index in range(2):
            for evt_index in range(len(self.event_energies[det_index])):
                erg = self.event_energies[det_index][evt_index]
                time = self.event_times[det_index][evt_index] - self.beam_center_times[det_index]
                eff = eff_lists[det_index][evt_index]
                dt_corr = self.dead_time_corrs[det_index][evt_index]

                br_energy[0] = erg
                br_time[0] = time
                br_eff[0] = eff
                br_dead_corr[0] = dt_corr

                self.trees[det_index].Fill()

                prog_index += 1
                prog.log(prog_index)

        [tree.Write() for tree in self.trees]

    def __init__(self, list_file_maestro: MaestroListFile, list_file_mpant: MPANTList,
                 eff_fit_maestro: ModelResult, eff_fit_mpant: ModelResult, new_name):
        """

        Args:
            list_file_maestro: Listfile class from LLNL detector
            list_file_mpant: Listfile class from IAC detector
            eff_fit_maestro: Corresponding efficiency for LLNL detector
            eff_fit_mpant: See eff_fit_maestro
        """

        root_path, pickle_path, marshal_path = self.__get_paths__(new_name)

        if not root_path.parent.exists():
            root_path.parent.mkdir()

        assert isinstance(eff_fit_maestro, ModelResult)
        assert isinstance(eff_fit_mpant, ModelResult)
        self.eff_fits: List[ModelResult] = [eff_fit_maestro, eff_fit_mpant]

        self.beam_center_times, self.beam_durations, self.beam_abs_times =\
            IACShot.get_beam_times(list_file_maestro, list_file_mpant)
        self.trees = None

        self.erg_centers = [list_file_maestro.erg_centers, list_file_mpant.erg_centers]
        self.erg_bins = [list_file_maestro.erg_bins, list_file_mpant.erg_bins]
        self.event_times = [list_file_maestro.times, list_file_mpant.times]
        self.event_energies = [list_file_maestro.energies, list_file_mpant.energies]
        self.dead_time_corrs = [list_file_maestro.deadtime_corrs, list_file_mpant.deadtime_corrs]

        self.pickle(pickle_path, marshal_path)
        root_file = ROOT.TFile(str(root_path), 'recreate')
        IACShot.root_files.append(root_file)
        self._build_tree()

    @cached_property
    def integrated_spectra(self):
        return [np.array([len(evts) for evts in det_list]) for det_list in self._list_of_data]

    @cached_property
    def _list_of_data(self):
        _list_of_data = [[], []]

        # append MT list corresponding to each channel (for both detectors)
        [_list_of_data[0].append([]) for _ in range(len(self.erg_centers[0]))]
        [_list_of_data[1].append([]) for _ in range(len(self.erg_centers[1]))]

        for det_index in range(2):
            for erg, time in zip(self.event_energies[det_index], self.event_times[det_index]):
                channel_index = self.erg_bin_index(erg, det_index)
                _list_of_data[det_index][channel_index].append(time)
        return _list_of_data

    def pickle(self, pickle_path, marshal_path):
        with open(pickle_path, 'wb') as f:  # pickle things that cannot be marshal'd
            pickle_data = {'eff_fits': self.eff_fits}
            dill.dump(pickle_data, f)

        with open(marshal_path, 'wb') as f:
            d = {'event_energies': [np.array(l, dtype=np.float64) for l in self.event_energies],
                 'event_times': [np.array(l, dtype=np.float64) for l in self.event_times],
                 'beam_abs_times': [np.array(l, dtype=np.float64) for l in self.beam_abs_times],
                 'erg_centers': [np.array(l, dtype=np.float64) for l in self.erg_centers],
                 'erg_bins': [np.array(l, dtype=np.float64) for l in self.erg_bins]}

            marshal.dump(d, f)

    @classmethod
    def load(cls, name) -> IACShot:
        self = cls.__new__(cls)
        root_path, pickle_path, marshal_path = self.__get_paths__(name)

        with open(pickle_path, 'rb') as f:
            pickle_data = dill.load(f)
            self.eff_fits = pickle_data['eff_fits']

        with open(marshal_path, 'rb') as f:
            marshal_data = marshal.load(f)
            for k, v in marshal_data.items():
                setattr(self, k, [np.frombuffer(b, dtype=np.float64) for b in v])

        root_file = ROOT.TFile(str(root_path))
        IACShot.root_files.append(root_file)
        self.trees = root_file.Get('tree1'), root_file.Get('tree2')
        return self

    def __repr__(self):
        assert False

    def erg_bin_index(self, erg, det_index):
        return np.searchsorted(self.erg_bins[det_index], erg, side='right') - 1

    def eval_efficiency(self, erg, det_index, return_ufloat=False):
        fit = self.eff_fits[det_index]

        if not hasattr(erg, '__iter__'):
            out = fit.eval(x=[erg])
            err = fit.eval_uncertainty(x=[erg])
            if return_ufloat:
                return ufloat(out[0], err[0])
        else:
            out = fit.eval(x=erg)
            err = fit.eval_uncertainty(x=erg)
            if return_ufloat:
                return unp.uarray(out, err)

    def plot_integrated_spectra(self, eff_correct=True, use_one_ax=None, erg_range=None):
        if erg_range is None:
            erg_range = min([self.erg_bins[0][0], self.erg_bins[1][0]]),\
                        max([self.erg_bins[0][-1], self.erg_bins[1][-1]])
        if use_one_ax is None:
            fig, axs = plt.subplots(2, 1, sharex="all", figsize=(9, 9*9/16))
            axs = axs.flatten()
        else:
            if not isinstance(use_one_ax, type(plt.gca())):
                plt.figure()
                use_one_ax = plt
            axs = [use_one_ax]*2
        labels = ['Upstream', 'downstream']
        for i, ax in enumerate(axs):
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            result = unp.uarray(self.integrated_spectra[i], np.sqrt(self.integrated_spectra[i]))
            min_i, max_i = np.searchsorted(self.erg_centers[i], erg_range, side='right')
            if min_i == max_i:
                continue
            result = result[min_i: max_i]
            if eff_correct:
                eff = self.eval_efficiency(self.erg_centers[i][min_i: max_i], i, True)
                result = result/eff
            y = unp.nominal_values(result)
            yerr = unp.std_devs(result)
            mpl_hist(self.erg_bins[i][min_i: max_i+1], y, yerr, ax=ax, label=labels[i])
            ax.legend()

        if not use_one_ax:
            plt.subplots_adjust(hspace=0.16)
            return axs
        else:
            return use_one_ax

    def time_dependence_over_range(self, min_erg, max_erg, det_index, time_bin_edges, eff_correct=True,
                                  ):
        """
        Returns eff. corrected time dependence in a energy range.
        Args:
            erg: Todo
            det_index: Upstream or downstream detector? (0 or 1)
            bins: Just like np.histogram
            account4eff: Correct for efficiency?
            return_ufloat: Calc. error?

        Returns:

        """
        time_bin_edges = np.array(time_bin_edges)
        bin_widths = time_bin_edges[1:] - time_bin_edges[:-1]
        min_channel = self.erg_bin_index(min_erg, det_index)
        max_channel = self.erg_bin_index(max_erg, det_index)
        weights = np.zeros(len(time_bin_edges)-1)
        sum_w2 = np.zeros_like(weights)
        weights_errors = np.zeros_like(weights)
        for ch_index in range(min_channel, min([max_channel+1, len(self._list_of_data[det_index])])):
            times = self._list_of_data[det_index][ch_index]
            times = np.array(list(filter(lambda x: time_bin_edges[0] <= x < time_bin_edges[-1], times)))
            if not len(times):
                continue
            time_is = np.searchsorted(time_bin_edges, times, side='right') - 1
            if eff_correct:
                eff = self.eval_efficiency(self.erg_centers[det_index][ch_index], det_index, True)
                weights[time_is] += 1/eff.n
                sum_w2[time_is] += 1/eff.n**2
                weights_errors[time_is] += eff.n
            else:
                weights[time_is] += 1
        weights = unp.uarray(weights, np.sqrt(sum_w2))
        weights += unp.uarray(np.zeros_like(weights_errors), weights_errors)
        print(bin_widths)
        weights /= bin_widths  # counts -> Hz

        return weights

    def plot_time_dependence(self, min_erg=None, max_erg=None, num_time_bins=None, eff_correct=True):
        """

        Args:
            min_erg:
            max_erg:
            num_time_bins: If None, auto.
            eff_correct: Whether or not to perform efficiency correction.

        Returns:

        """
        if min_erg is None:
            min_erg = min([self.erg_centers[0][0], self.erg_centers[1][0]])
        if max_erg is None:
            max_erg = max([self.erg_centers[0][-1], self.erg_centers[1][-1]])
        assert max_erg >= min_erg
        fig, axs = plt.subplots(2, 2, figsize=(9, 9 * 9 / 16))
        gs = axs[-1][-1].get_gridspec()
        for ax in axs[:, -1]:
            ax.remove()
        erg_axs = fig.add_subplot(gs[:, -1])
        # print(axs[1,0].remove())
        time_axs = axs[:, 0]

        labels = ['Upstream', 'downstream']
        erg_axs.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        erg_window_width = max_erg-min_erg
        ax = self.plot_integrated_spectra(use_one_ax=erg_axs, erg_range=[min_erg-erg_window_width, max_erg+erg_window_width],
                                          eff_correct=eff_correct)

        ax.fill_between([min_erg, max_erg], [ax.get_ylim()[0]] * 2, [ax.get_ylim()[1]] * 2, alpha=0.4, color='red',
                        label='energy window')
        ax.legend()
        if num_time_bins is None:
            time_bins = np.histogram_bin_edges(np.concatenate(self.event_times))
        else:
            raise NotImplementedError("Make more better")

        for i, ax in enumerate(time_axs):
            y = self.time_dependence_over_range(min_erg, max_erg, i, time_bins,
                                                eff_correct=eff_correct)
            yerr = unp.std_devs(y)
            y = unp.nominal_values(y)
            mpl_hist(time_bins, y, yerr, ax=ax, label=labels[i])
            ax.set_xlabel("Time since beam fire [s]")
            ax.set_ylabel("Counts/s [Hz]")



def load_efficiency(rel_data_dir_path) -> ModelResult:
    """
    Load efficiency from cal_data/rel_data_dir_path.
    Saved efficiencies are always named "cal.pickle", so just get `rel_data_dir_path` right.
    See eff_cal.py for creating efficiency fits.

    Args:
        rel_data_dir_path:

    Returns:

    """
    path = top_level_data_path/rel_data_dir_path/'cal.pickle'
    with open(path, 'rb') as f:
        fit: ModelResult = pickle.load(f)
    return fit


if __name__ == '__main__':
    #       todo:
    #        - Make a way to split Maestro and MPANT list files for multiple shots in short succession
    #        -
    #
    #
    #        interpolate percent live time.
    # with open('delete.marshal', 'wb') as f:
    #     d = {'a': [[1,2,3,4,5], [1,2.]], 'b':[[],[],[],[],[],[]]}
    #     marshal.dump(d, f)
    # with open('delete.marshal', 'rb') as f:
    #     print(marshal.load(f))
    # from JSB_tools import Nuclide

    list_file_mpant = MPANTList("/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/list_files/jeff002.txt",primary_adc_number=1, aux_adc_number=3)
    list_file_maestro = MaestroListFile("/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/list_files/test_beam.Lis")
    # list_file_maestro.plot_count_rate()
    shot = IACShot(list_file_maestro, list_file_mpant, load_efficiency('our_det/2021-08-17'), load_efficiency('our_det/2021-08-17'),
                   'test')
    # shot = IACShot.load('test')
    # plt.figure()
    # shot.plot_integrated_spectra(eff_correct=False, use_one_ax=plt.gca(), erg_range=[1000, 1500])
    # with plt.xkcd():
    #     shot.plot_time_dependence()
    # plt.show()
    tb = ROOT.TBrowser()

    ROOT_loop()


