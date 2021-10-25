import matplotlib.pyplot as plt

from JSB_tools.list_reader import MaestroListFile
from JSB_tools.spe_reader import SPEFile
import timeit
import cProfile
from lmfit.models import GaussianModel, LinearModel
from lmfit.model import CompositeModel
import numpy as np
import pickle

bins = np.linspace(-10 + 3, 10 + 3)
x = (bins[1:] + bins[:-1])/2
points = 2*np.random.randn(10000) + 3
y_gauss, _ = np.histogram(points, bins=bins)
yerr_gauss = np.sqrt(y_gauss)
weights = np.where(yerr_gauss > 0, yerr_gauss, 1)
weights = 1.0/weights

y_lin = (2*x + 5) * np.random.uniform(0.9, 1.1, len(x))

model_lin = LinearModel()
fit_result_lin = model_lin.fit(data=y_lin, x=x)
fit_result_lin.plot_fit()


model = GaussianModel()
params = model.guess(x=x, data=y_gauss, weights=weights)
fit_result_guass = model.fit(x=x, data=y_gauss, weights=weights)

# plt.errorbar(x, y_gauss, yerr=yerr_gauss)
with open("test.pickle", 'wb') as f:
    pickle.dump(fit_result_guass, f)
    print(type(fit_result_guass))
#
# with open("test.pickle", 'rb') as f:
#     fit_result = pickle.load(f)
# plt.figure()
# fit_result.plot_fit()


plt.show()




# spe = MaestroListFile('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/shot133.Lis')
# spe.build_spe()
# cProfile.run('spe.__build_spe__()', )
# t1 = timeit.timeit('spe.__build_spe__()', number=100, globals=globals())
# print(t1)

