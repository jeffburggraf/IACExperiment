import numpy as np
from matplotlib import pyplot as plt
from JSB_tools import mpl_hist, Nuclide
from pathlib import Path
from analysis import Shot
from JSB_tools import convolve_gauss
from lmfit.models import GaussianModel, ExponentialModel
from uncertainties import unumpy as unp, ufloat


# ========================================================
nuclides = ['Sb132']   # ['Sr94', 'Xe139']
fit_width = 2
gamma_indexes = {}
fix_half_life = ['Sr93']
# ========================================================
yields = {}


def get_fit_spec(e, w):
    b0 = np.searchsorted(bins,e - w/2, side='right') - 1
    b1 = np.searchsorted(bins, e + w/2, side='right') - 1
    return x_all[b0:b1], y_all[b0:b1], b0, b1


shots = Shot.find_shots(flow='100001',
                        tube_len=6.14,
                        flow_stop=None,
                        foil_pos='upstream',
                        beam_duration=3,
                        mylar=0,
                        cold_filter=False,
                        eval_func="self.he_flow == self.ar_flow and self.ar_flow in [0.5, 0.25]")
for n in nuclides:
    gamma_index = gamma_indexes.get(n, 0)
    n = Nuclide.from_symbol(n)
    gamma_line = n.decay_gamma_lines[gamma_index]
    gamma_erg = gamma_line.erg.n
    peak_width =  Shot.get_peak_width(gamma_erg)
    print(n)

    for shot in shots:
        fig, axs = plt.subplots(1, 2)
        ax1, ax2 = axs
        list_file = shot.list
        y_all, bins = list_file.get_erg_spectrum(erg_min=80, eff_corr=True, return_bin_edges=True, remove_baseline=True)

        x_all = 0.5*(bins[1:] + bins[:-1])
        b_widths = bins[1:] - bins[:-1]

        y_all /= b_widths

        x_fit, y_fit, b0, b1 = get_fit_spec(gamma_erg, peak_width * 1.1)

        y_plot = y_all[b0-20: b1+20]
        x_plot = x_all[b0-20: b1+20]

        weights = 1.0/np.where(unp.std_devs(y_fit) > 0, unp.std_devs(y_fit), 1)

        model = GaussianModel()

        params = model.guess(unp.nominal_values(y_fit), x=x_fit)
        params['center'].value = gamma_line.erg.n
        params['sigma'].value = Shot.get_sigma(gamma_line.erg.n)

        fit_result = model.fit(data=unp.nominal_values(y_fit), x=x_fit, params=params, weights=weights)

        ax1.plot(x_fit, fit_result.eval(params=params), label='model')
        ax1.errorbar(x_plot, unp.nominal_values(y_plot), unp.std_devs(y_plot), label='data')

        # time_b_width = n.half_life.n / 5
        # time_bins = np.linspace(0, shot.max_time, int(shot.max_time // time_b_width))
        n_counts = np.sum(list_file.get_erg_spectrum(gamma_erg - peak_width / 2, gamma_erg + peak_width / 2,
                                                     nominal_values=True, eff_corr=False, remove_baseline=True))
        time_bin_width = min([n.half_life.n / 12, 15])
        time_bin_width = max([3, time_bin_width])
        time_bins = np.arange(0, shot.max_time + time_bin_width, time_bin_width)

        time_dep, _, _ = list_file.get_time_dependence(gamma_erg,
                                                       signal_window_kev=Shot.get_peak_width(gamma_erg),
                                                       eff_corr=True,
                                                       make_rate=True,
                                                       bins=time_bins,
                                                       nominal_values=False)
        temp_y = convolve_gauss(unp.nominal_values(time_dep), 2)
        max_i = np.argmax(temp_y)
        fit_times = 0.5*(time_bins[1:] + time_bins[:-1])
        time_of_max_Rate = fit_times[max_i]

        fit_times = fit_times[max_i:]
        fit_times -= fit_times[0]
        fit_rates = time_dep[max_i:]
        fit_rates_weights = 1.0/np.where(unp.std_devs(fit_rates)>0, unp.std_devs(fit_rates), 1)
        fit_rates = unp.nominal_values(fit_rates)
        model = ExponentialModel()
        params = model.guess(data=fit_rates, x=fit_times)
        # print(params)
        params['decay'].set(value=1.0/n.decay_rate.n)

        if n.name in fix_half_life:
            params['decay'].set(min=0.8*params['decay'].value, max=1.2*params['decay'].value)

        print(params)

        fit_result = model.fit(x=fit_times, data=fit_rates, weights=fit_rates_weights, params=params)
        params = fit_result.params
        # print(fit_result.fit_report())
        print(f"t_1/2 fit: {0.693*fit_result.params['decay'].value}")
        print("After acq. correction: ", 1.0/(1-np.e**(-n.decay_rate.n*list_file.times[-1])))
        fit_plot_x = np.linspace(0, time_bins[-1], 200)
        fit_plot_y = fit_result.eval(params=params, x=fit_plot_x-time_of_max_Rate)
        fit_plot_y_err = fit_result.eval_uncertainty(params=params, x=fit_plot_x-time_of_max_Rate)
        ax2.plot(fit_plot_x, fit_plot_y, label='Fit')

        mpl_hist(time_bins, time_dep, ax=ax2)

        n_decays = ufloat(fit_result.params['amplitude'].value, fit_result.params['amplitude'].stderr)
        n_decays /= gamma_line.intensity

        fig.suptitle(f"{n.name},   $t_{{1/2}}$={n.half_life.n}\n"
                      f"Shot {shot.shotnum}; "
                      f"tmax: {shot.max_time:.1f} or {shot.max_time / n.half_life.n:.1f}$t_{{1/2}}$;"
                      f" n_decays:{n_decays:.2e}")

        attribs = {'cold': shot.cold_filter, 'flow': [shot.he_flow, shot.ar_flow], }

        print(f"\tyield: {n_decays}, shot{shot.shotnum}, {attribs}")

plt.show()








