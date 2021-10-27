import matplotlib.pyplot as plt

from JSB_tools.list_reader import MaestroListFile
from JSB_tools.spe_reader import SPEFile
import timeit
import cProfile
from lmfit.models import GaussianModel, LinearModel
from lmfit.model import CompositeModel
import numpy as np
import pickle


class Mixin:
    def __init__(self):
        pass

    def f(self):
        print(self.path)


class A(Mixin):
    def __init__(self, path):
        super().__init__()
        self.path = path


a = A('omg')
a.f()







# spe = MaestroListFile('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/shot133.Lis')
# spe.build_spe()
# cProfile.run('spe.__build_spe__()', )
# t1 = timeit.timeit('spe.__build_spe__()', number=100, globals=globals())
# print(t1)

