from pathlib import Path
import ROOT
import numpy as np
from JSB_tools.list_reader import MaestroListFile
from matplotlib import pyplot as plt
from typing import Union
from mpant_reader import MPANTList


list_file_dir = Path.cwd()/'exp_data'/'list_files'


def get_ortec_beam_time(l: Union[MaestroListFile, MPANTList]):
    if isinstance(l, MaestroListFile):
        beam_on_times = l.realtimes[np.where(l.sample_ready_state == 1)]
        beam_on_center_time = np.median(beam_on_times)
        beam_duration = 2 * np.percentile(beam_on_times - beam_on_center_time, 98)
        return beam_on_center_time, beam_duration
    else:
        raise NotImplementedError


class IACShot:
    save_path = Path.cwd()/'exp_data'/'packaged_shots'

    @classmethod
    def load(cls, name, new_name=None):
        pass

    def __init__(self, list_file: Union[MaestroListFile, MPANTList], efficiency_function, new_name=None):
        if new_name is None:
            new_name = list_file.file_name

        root_file = ROOT.TFile(str(IACShot.save_path), 'recreate')
        tree = ROOT.TTree('tree', 'tree')
        if isinstance(list_file, MaestroListFile):
            self.beam_center_time, self.beam_duration = get_ortec_beam_time(list_file)
        else:
            raise NotImplementedError

        br_time = np.array([0.])
        br_energy = np.array([0.])
        br_eff = np.array([0.])
        br_dead_corr = np.array([0.])

        tree.Branch("t", br_time, 't/F')
        tree.Branch("erg", br_energy, 'erg/F')
        tree.Branch("eff", br_eff, 'eff/F')
        tree.Branch('dead_corr', br_dead_corr, 'dead_corr/F')
        print(len(list_file.percent_live), len(list_file.times))

#          todo: interpolate percent live time
#           Add class for MPANT List file











if __name__ == '__main__':
    list_file = MaestroListFile(list_file_dir / 'test.Lis')
    shot = IACShot(list_file, lambda x:x)

    list_file.plot_sample_ready()
    plt.show()



