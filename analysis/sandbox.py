import warnings
from pathlib import Path
import re
from matplotlib import pyplot as plt
import numpy as np
from typing import Tuple
from JSB_tools.list_reader import MaestroListFile
from mpant_reader import MPA
from pickle_shots import get_maesto_list_shot_paths, get_mpant_mca_shot_paths

shot_dirs = []

mca_paths = get_mpant_mca_shot_paths()
maestro_list_paths = get_maesto_list_shot_paths()


class Shot:
    def __init__(self, shot_num):
        self.shot_num = int(shot_num)
        try:
            self.mca = MPA(mca_paths[self.shot_num])
        except (FileNotFoundError, KeyError):
            warnings.warn(f"No MPA data for shot {shot_num}")
            self.mca = None
        try:
            self.maestro_list = MaestroListFile.from_pickle(maestro_list_paths[shot_num])
        except (FileNotFoundError, KeyError):
            raise

        self.maestro_zero_time = np.median(self.maestro_list.realtimes[np.where(self.maestro_list.sample_ready_state == 1)])
        self.maestro_list.times = self.maestro_list.times - self.maestro_zero_time
        self.maestro_list.__needs_updating__ = True

    def get_binned_times(self, time_bins, erg_min, erg_max, maestro=True) -> Tuple[np.ndarray, np.ndarray]:
        """
        Returns an array of the number of counts in each time bin subject to an energy cut.
        Args:
            time_bins:
            erg_min:
            erg_max:
            maestro:

        Returns: Bin values,  bin centers

        """
        time_bins = np.array(time_bins)
        if maestro:
            times = self.maestro_list.get_times_in_range(erg_min=erg_min, erg_max=erg_max)
            values, _ = np.histogram(times, bins=time_bins)
            b_centers = (time_bins[1:] + time_bins[:-1])/2
            return b_centers, values
        else:
            raise NotImplementedError

    def plot_218(self, center=):



if __name__ == '__main__':
    s = Shot(3)
    t_bins = np.arange(0, 200, 3)
    x, y = s.get_binned_times(t_bins, 217, 220)
    x1, y1 = s.get_binned_times(t_bins, 214, 217)

    plt.plot(x, y)
    plt.plot(x1, y1, label='off-peak')
    plt.legend()

    s.maestro_list.plot_erg_spectrum()

    plt.show()