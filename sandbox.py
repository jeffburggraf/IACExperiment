import matplotlib.pyplot as plt

from JSB_tools.list_reader import MaestroListFile
from JSB_tools.spe_reader import SPEFile
import timeit
import cProfile
from lmfit.models import GaussianModel, LinearModel
from lmfit.model import CompositeModel
import numpy as np
import pickle

from analysis import Shot
from FFStopping import GetPoints
for shot in Shot.find_shots(foil_pos='upstream', beam_duration=3):
    print(shot.shotnum)







# spe = MaestroListFile('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/shot133.Lis')
# spe.build_spe()
# cProfile.run('spe.__build_spe__()', )
# t1 = timeit.timeit('spe.__build_spe__()', number=100, globals=globals())
# print(t1)

