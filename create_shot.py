from pathlib import Path
import ROOT
from JSB_tools import ROOT_loop
import numpy as np
from JSB_tools.list_reader import MaestroListFile
from matplotlib import pyplot as plt
from typing import Union
from mpant_reader import MPANTList


list_file_dir = Path.cwd()/'exp_data'/'list_files'


def get_beam_time(l: Union[MaestroListFile, MPANTList]):
    if isinstance(l, MaestroListFile):
        beam_on_times = l.realtimes[np.where(l.sample_ready_state == 1)]
    elif isinstance(l, MPANTList):
        beam_on_times = l.times(adc=2)
    else:
        assert False

    beam_on_center_time = np.median(beam_on_times)
    beam_duration = 2 * np.percentile(beam_on_times - beam_on_center_time, 98)
    return beam_on_center_time, beam_duration


root_files = []

class IACShot:
    save_path = Path.cwd()/'exp_data'/'packaged_shots'

    @classmethod
    def load(cls, name, new_name=None):
        pass

    def __init__(self, list_file: Union[MaestroListFile, MPANTList], efficiency_function, new_name=None):
        if new_name is None:
            new_name = list_file.file_name
        root_path = (IACShot.save_path/new_name).with_suffix(".root")
        if not root_path.parent.exists():
            root_path.parent.mkdir()

        root_file = ROOT.TFile(str(root_path), 'recreate')

        tree = ROOT.TTree('tree', 'tree')
        self.beam_center_time, self.beam_duration = get_beam_time(list_file)

        br_time = np.array([0.], dtype=np.float32)
        br_energy = np.array([0.], dtype=np.float32)
        br_eff = np.array([0.], dtype=np.float32)
        br_dead_corr = np.array([0.], dtype=np.float32)

        tree.Branch("t", br_time, 't/F')
        tree.Branch("erg", br_energy, 'erg/F')
        tree.Branch("eff", br_eff, 'eff/F')
        tree.Branch('dead_corr', br_dead_corr, 'dead_corr/F')

        if isinstance(list_file, MaestroListFile):
            event_times = list_file.times
            event_energies = list_file.energies
        elif isinstance(list_file, MPANTList):
            event_times = list_file.times(adc=1)
            event_energies = list_file.energies(adc=1)
        else:
            assert False

        # Todo: Pickle eff function, listfile, beam_center_time, whatever else.

        for i in range(len(event_energies)):
            erg = event_energies[i]
            time = event_times[i] - self.beam_center_time

            eff = efficiency_function(erg)
            br_energy[0] = erg
            br_time[0] = time
            br_eff[0] = eff
            tree.Fill()

        tree.Write()
        root_files.append(root_file)


        # root_file.Write()


#          todo: interpolate percent live time
#           Add class for MPANT List file











if __name__ == '__main__':

    list_file = MPANTList("/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/list_files/beamgun003.txt")

    shot = IACShot(list_file, lambda x:x)

    tb = ROOT.TBrowser()

    ROOT_loop()


