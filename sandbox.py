# from JSB_tools.list_reader import MaestroListFile
# from JSB_tools import mpl_hist
# from analysis import time_offset
# import matplotlib.pyplot as plt
# l = MaestroListFile.from_pickle('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/shot132.Lis')
# time_offset(l)
#
# ys, yb, bins = l.get_time_dependence(1427.7, bins='auto')
#
# ax = mpl_hist(bins, ys)
# mpl_hist(bins, yb, ax=ax)
# plt.show()

instnaces = {}




@cached_cls
class A:
    def __init__(self, a):
        self.a = a


a = A(3)
a = A(2)
a = A(2)
