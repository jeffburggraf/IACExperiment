from pathlib import Path
import ROOT
from JSB_tools import ROOT_loop, TBrowser
import numpy as np
from JSB_tools.list_reader import MaestroListFile
from matplotlib import pyplot as plt
from typing import Union
from mpant_reader import MPANTList
import dill
import marshal


list_file_dir = Path.cwd()/'exp_data'/'list_files'


def get_beam_time(l: Union[MaestroListFile, MPANTList]):
    if isinstance(l, MaestroListFile):
        beam_abs_times = l.realtimes[np.where(l.sample_ready_state == 1)]
    elif isinstance(l, MPANTList):
        beam_abs_times = l.get_times(adc=2)
    else:
        assert False

    beam_on_center_time = np.median(beam_abs_times)
    beam_duration = 2 * np.percentile(beam_abs_times - beam_on_center_time, 98)
    return beam_on_center_time, beam_duration, beam_abs_times


class IACShot:
    """
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
        out.eff_function = lambda x: 1
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

        self.tree = ROOT.TTree('tree', 'tree')
        self.tree.Branch("t", br_time, 't/F')
        self.tree.Branch("erg", br_energy, 'erg/F')
        self.tree.Branch("eff", br_eff, 'eff/F')
        self.tree.Branch('dead_corr', br_dead_corr, 'dead_corr/F')

        for i in range(len(self.event_energies)):
            erg = self.event_energies[i]
            time = self.event_times[i] - self.beam_center_time

            eff = self.eff_function(erg)
            br_energy[0] = erg
            br_time[0] = time
            br_eff[0] = eff
            self.tree.Fill()

        self.tree.Write()

    def __init__(self, list_file: Union[MaestroListFile, MPANTList], eff_function, new_name=None):

        if new_name is None:
            new_name = list_file.file_name

        root_path, pickle_path, marshal_path = self.__get_paths__(new_name)

        if not root_path.parent.exists():
            root_path.parent.mkdir()

        assert hasattr(eff_function, '__call__')
        self.eff_function = eff_function

        self.beam_center_time, self.beam_duration, self.beam_abs_times = get_beam_time(list_file)
        self.tree = None

        self._list_of_data = []

        self.event_times = list_file.times
        self.event_energies = list_file.energies

        [self._list_of_data.append([]) for i in range(len(list_file.erg_bins) - 1)]  # append MT list 4 each erg ch
        for erg, time in zip(self.event_energies, self.event_times):  # fill self._list_of_data
            i = list_file.erg_bin_index(erg)
            self._list_of_data[i].append(time)

        self.pickle(pickle_path, marshal_path)
        root_file = ROOT.TFile(str(root_path), 'recreate')
        IACShot.root_files.append(root_file)
        self._build_tree()

    def pickle(self, pickle_path, marshal_path):
        with open(pickle_path, 'wb') as f:  # pickle things that are not computationally intensive
            pickle_data = {'beam_abs_times': self.beam_abs_times, 'eff_function': self.eff_function}
            dill.dump(pickle_data, f)

        with open(marshal_path, 'wb') as f:
            marshal.dump(self._list_of_data, f)

    @classmethod
    def load(cls, name):
        self = cls.__new__(cls)
        root_path, pickle_path, marshal_path = self.__get_paths__(name)

        with open(pickle_path, 'rb') as f:
            pickle_data = dill.load(f)
            self.beam_abs_times = pickle_data['beam_abs_times']
            self.eff_function = pickle_data['eff_function']

        with open(marshal_path, 'rb') as f:
            self._list_of_data = marshal.load(f)

        root_file = ROOT.TFile(str(root_path))
        IACShot.root_files.append(root_file)
        self.tree = root_file.Get('tree')


if __name__ == '__main__':
    #       todo: interpolate percent live time.
    from JSB_tools import Nuclide
    n = Nuclide.from_symbol('Na22')
    for g in n.decay_gamma_lines:
        print(g)
    # shot = IACShot.gen_fake_data()
    # shot = IACShot.load('fake')
    # list_file = MPANTList("/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/list_files/beamgun003.txt")
    # list_file2 = MaestroListFile("/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/list_files/BeamGunSimulation.Lis")
    # shot = IACShot(list_file2, lambda x:x)

    #
    # tb = ROOT.TBrowser()

    # ROOT_loop()


