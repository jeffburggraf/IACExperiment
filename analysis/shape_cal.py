import pickle
import matplotlib.pyplot as plt
import numpy as np
from lmfit.models import GaussianModel
from analysis import Shot
from JSB_tools import mpl_hist
from scipy.signal import find_peaks
from uncertainties import unumpy as unp, ufloat
from JSB_tools.spe_reader import SPEFile
from uncertainties.umath import sqrt as usqrt

f1 = SPEFile("/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/EffCalFriday/Eu152-0.Spe")
f1_ergs = [121.8, 244.7, 344.3, 779, 964,1112, 1408]
bins_1 = f1.erg_bins


f2 = SPEFile('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/EffCalFriday/Y88-0.Spe')
f2_ergs = [898.2, 1836, ]



def eval_fit(_x, coeffs):
    _x = np.array(_x)
    return np.sum([c*_x**i for i, c in enumerate(reversed(coeffs))], axis=0)


def get_approx_fwhm(erg):
    return eval_fit([erg], [-2.38396079e-08,  2.52368849e-04, 4.43790741e-01])[0]*2.35


def std(xs, weights):
    if not isinstance(xs, np.ndarray):
        xs = np.array(xs)

    if not isinstance(weights, np.ndarray):
        weights = np.array(weights)

    mean = np.sum(xs*weights)/np.sum(weights)
    var = np.sum(weights*(xs - mean)**2)/np.sum(weights)
    return usqrt(var)


ergs = []
sigmas = []

for peak_centers, f in [(f1_ergs, f1), (f2_ergs, f2)]:
    counts = f.counts
    assert isinstance(f, SPEFile)
    ax = f.plot_erg_spectrum()
    for peak in peak_centers:
        fit_width = get_approx_fwhm(peak)*2
        ergs.append(peak)
        model = GaussianModel()
        ilow = f.__erg_index__(peak - fit_width/2)
        ihigh = f.__erg_index__(peak + fit_width/2)
        fit_y = counts[ilow: ihigh]

        weights = 1.0/np.where(unp.std_devs(fit_y) > 0, unp.std_devs(fit_y), 1)
        fit_y = unp.nominal_values(fit_y)
        bins_fit = f.erg_bins_cut(peak - fit_width/2, peak + fit_width/2)
        fit_x = 0.5*(bins_fit[1:] + bins_fit[:-1])
        params = model.guess(x=fit_x, data=fit_y)
        fit_result = model.fit(params=params, x=fit_x, data=fit_y, weights=weights)
        eval_x = np.linspace(fit_x[0], fit_x[-1], 100)
        ax.plot(eval_x, fit_result.eval(params=params, x=eval_x))
        sigma = ufloat(fit_result.params['sigma'].value, fit_result.params['sigma'].stderr)
        sigmas.append(sigma)
        fit_center = params['center']
        ax.axvspan(fit_center-Shot.get_peak_width(peak)/2, fit_center+Shot.get_peak_width(peak)/2, alpha=0.25, color='black')

plt.figure()
fit_coefs = np.polyfit(ergs, unp.nominal_values(sigmas), deg=2, w=1.0/unp.std_devs(sigmas))
sigma_fit_ergs = np.linspace(0, 2000, 300)

plt.errorbar(ergs, unp.nominal_values(sigmas), unp.std_devs(sigmas), ls='None', marker='o', label='Fitted sigmas')
plt.plot(sigma_fit_ergs, eval_fit(sigma_fit_ergs, fit_coefs), label='Poly regression')
plt.ylabel('Sigma [KeV]')
plt.xlabel('Erg [KeV]')
plt.legend()

print('Shape cal. quadratic coeffs (erg -> sigma): ', list(reversed(fit_coefs)))
plt.show()


