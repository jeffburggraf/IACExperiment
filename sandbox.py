# import numpy as np
# from JSB_tools.spe_reader import SPEFile
# from matplotlib import pyplot as plt
from openpyxl import load_workbook

wb = load_workbook(filename = r'/Users/burggraf1/Downloads/IAC Run Spreadsheet (4).xlsx')
# sheet_ranges = wb[]

sheet = wb['Sheet1']
for x in (sheet['A']):
    print(x.value)


#
# from JSB_tools import Nuclide
# #
# n = Nuclide.from_symbol('Y88')
# print(n.decay_daughters)
# for g in n.decay_gamma_lines:
#     print(g)
# # # s = SPEFile('/Users/burggraf1/PycharmProjects/IACExperiment/cal_data/our_det/2021-08-17/Eu-152-0.Spe')
# #
# # s.get_spectrum_hist().plot()
# # plt.show()
# import pickle
#
# with open("/Users/burggraf1/PycharmProjects/IACExperiment/cal_data/our_det/2021-08-17/cal.pickle", 'rb') as f:
#     d = pickle.load(f)
# print(d.eval(x=[2]))