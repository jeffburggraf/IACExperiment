import re
import warnings

import matplotlib.pyplot as plt
from JSB_tools import DecayNuclide
from JSB_tools import Nuclide, FissionYields
from pathlib import Path
from JSB_tools.MCNP_helper.outp_reader import OutP
import scipy
from JSB_tools import mpl_style
from scipy import integrate
import numpy as np
import pickle
from uncertainties import unumpy as unp
"""
Interesting examples:
Xe140: Does not suffer from decay chain contributions. Nearly 100% direct from fission. 
Rb89 - rising to its peak very late in time. 
La-144 - 397 KeV
"""
c_per_second = (192/3.0)*1E-6
charge_per_electron = 1.602E-19
n_electrons = 3*c_per_second / charge_per_electron

# =======================================================
# List of ffs (as strings), or an int to plot top N fission fragments.
ffs = 'Nb102'
# mpl_style(usetex=True, fontscale=1.)
save_fig = False
decay_rate = True  # If True, plot decay rate instead of nuclei number.
rel_yield_thresh = 0.01  # yields with integrated fractions of total below this value are not plotted.
arb_yaxis = False  # normalizes the y to smaller numbers.
fraction_escape = 0.1
t_max = 600
t_b_width = 1
plot_integral_yield = False
figsize = (8, 8)  # default = (12, 8)
fast = False  # for debugging--runs faster
# =======================================================

if isinstance(ffs, str):
    ffs = [ffs]


if hasattr(ffs, '__iter__'):
    if not len(ffs) == 1:
        save_fig = False

outp = OutP(Path(__file__).parent.parent/'mcnp'/'sims'/'du_shot131'/'outp')
tally_down = outp.get_f4_tally('Active down')
u238 = Nuclide.from_symbol('U238')

yields = FissionYields('U238', 'gamma', tally_down.energies)

_weights = tally_down.cell.atom_density * tally_down.dx_per_src * n_electrons * \
           u238.gamma_induced_fiss_xs.interp(tally_down.energies) * fraction_escape

yields.weight_by_erg(_weights)

times = np.arange(0, t_max, t_b_width)

#  decay_rates has the form:
# {'nuclide1':
#   {'fission': [y1, y2, ..., yn],
#    'parent_nuclide1': [y1, y2, ..., yn],
#    'parent_nuclide2': [y1, y2, ..., yn]
#    },
#  'nuclide2': {...}
#  }
# Where keys are nuclides that can be detected when they decay, and the values are dictionaries that list the
# contributions to the key via the decay of each parent.
decay_rates = {}

max_yield = None
i = 0

if fast:
    _yields_loop = list(yields.yields.items())[:20]
else:
    _yields_loop = yields.yields.items()

for parent, yield_ in _yields_loop:

    yield_ = sum(unp.nominal_values(yield_))

    if max_yield is None:
        max_yield = yield_
    elif yield_ < max_yield*0.001:
        break

    f = DecayNuclide(parent)
    decay_yields = f(times, decay_rate=decay_rate)

    if parent not in decay_rates:
        decay_rates[parent] = {}

    decay_rates[parent]['fission'] = yield_ * decay_yields[parent]
    del decay_yields[parent]

    for daughter, v in decay_yields.items():
        if daughter not in decay_rates:
            decay_rates[daughter] = {}

        try:
            decay_rates[daughter][parent] += yield_ * v
        except KeyError:
            decay_rates[daughter][parent] = yield_ * v
    i += 1


def srt(dict_):
    out = {'fission': dict_['fission']}
    del dict_['fission']

    for k, v in sorted(dict_.items(), key=lambda k_v: -max(k_v[1])):
        out[k] = v
    return out


if __name__ == '__main__':
    yields.undo_weighting()
    from JSB_tools import TabPlot

    tot_decay_yields = {}
    for k, v in decay_rates.items():
        tot = np.sum(list(v.values()), axis=0)
        tot_decay_yields[k] = np.trapz(tot, times)

    # sort all FFs by total decay yield
    tot_decay_yields = {k: v for k, v in sorted(tot_decay_yields.items(), key=lambda x: -x[-1])}

    if isinstance(ffs, int):  # select top N FFs
        ffs = [k for k in list(tot_decay_yields.keys())[:ffs]]

    yscale = None  # for if arb units are used
    tot_yield_scale = None

    if not save_fig:
        tab_plot = TabPlot(figsize=figsize)
    else:
        tab_plot = None

    gamma_printed = []

    def gamma_print(_n):
        if not isinstance(_n, Nuclide):
            _n = Nuclide.from_symbol(_n)
        if _n.name in gamma_printed:
            return
        print(f"{_n}:")
        for g in [g for g in _n.decay_gamma_lines if g.intensity.n > 0.01]:
            print(f'\t{g}')

        gamma_printed.append(_n.name)

    for ff in ffs:
        ff = Nuclide.from_symbol(ff)
        if not save_fig:
            ax = tab_plot.new_ax(ff.name)
        else:
            fig, ax = plt.subplots()

        if arb_yaxis:
            tot_yield_scale = 100/list(tot_decay_yields.values())[0]
        else:
            tot_yield_scale = 1

        try:
            tot_indep_yield = np.trapz(decay_rates[ff.name]['fission'], times)

        except KeyError:
            warnings.warn(f"No fission yield for {ff.name}... skipping")
            continue
            # assert False

        gamma_print(ff)

        if not plt.rcParams['text.usetex']:
            ax.set_title(f"{ff.name} - {ff.half_life:.1f} s\n"
                         f"Raw indep. yield: {np.average(yields[ff.name], weights=_weights)}")

        sorted_decay_sources = srt(decay_rates[ff.name])  # parent nuclei sosrted according to rates.
        ys = np.array(list(sorted_decay_sources.values()))
        cum_rates = np.cumsum(ys, axis=0)
        max_plt_value = max(cum_rates[-1])
        y1 = np.zeros_like(times)

        if yscale is None:
            if arb_yaxis:
                yscale = 100/max_plt_value
            else:
                yscale = 1

        cum_rates *= yscale
        ys *= yscale
        max_plt_value *= yscale

        handles, labels = [], []

        direct_handle, direct_label = [], []

        tot_decay_yield = np.trapz(cum_rates[-1], times)

        for index, parent_name in enumerate(sorted_decay_sources.keys()):
            y0 = y1
            y1 = cum_rates[index]
            dy = y1 - y0

            rel_yield = np.trapz(dy, times) / tot_decay_yield

            lbl = int(100 * rel_yield)
            if lbl == 0:
                lbl = f"{100*rel_yield:.1f}"
            extra_leg_label = f" ({lbl})"

            if parent_name != 'fission':
                legend_name = Nuclide.from_symbol(parent_name).latex_name
                linewidth = 1
                gamma_print(parent_name)
            else:
                legend_name = f'Direct-from-fission contribution'
                linewidth = 3

            if rel_yield < rel_yield_thresh and parent_name != 'fission':
                continue

            fill_leg_handle = ax.fill_between(times, y0, y1)
            line_leg_handle = ax.plot(times, y1, c='black', linewidth=linewidth)[0]

            label = f'{legend_name}{extra_leg_label}'

            if parent_name != 'fission':
                handles.append((fill_leg_handle, line_leg_handle))
                labels.append(label)
            else:
                direct_handle.append((fill_leg_handle, line_leg_handle))
                direct_label.append(label)

        if plot_integral_yield:
            ax_twin = ax.twinx()
            ax_twin.set_ylabel("Cumulative percent of tot. decays")
            if not save_fig:
                tab_plot.add_aux_axis(ax_twin)

            cum_yield = integrate.cumulative_trapezoid(cum_rates[-1], times, initial=0)
            cum_yield *= 1/max(cum_yield)
            ax_twin.plot(times, cum_yield, ls='--', label='Integral yield', c='black', alpha=0.7)
            ax_twin.legend(loc='center right')

        ax.set_xlabel("Time since pulse [s]")

        if decay_rate:
            ax.set_ylabel(f"Decay rate of {ff.latex_name} nuclei [arb. units]")
        else:
            ax.set_ylabel(f"# of {ff.name} nuclei [arb. units]")

        max_x = times[-np.searchsorted(cum_rates[-1][::-1], max(max_plt_value/150, cum_rates[-1][-2]))]

        ax.set_xlim(left=0, right=max_x)

        leg_title = r"Indirect contributions (\% of total)"

        leg = ax.legend(reversed(handles), reversed(labels), title=leg_title)
        ax.legend(direct_handle, direct_label, loc='center right')

        ax.add_artist(leg)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        if (not save_fig) and plot_integral_yield:
            all_axes = tab_plot.fig.get_axes()
            for axis in all_axes:
                legend = axis.get_legend()
                if legend is not None:
                    legend.remove()
                    all_axes[-1].add_artist(legend)

        if save_fig:
            path = f'transmutated_yields{ff.name}.png'
            if 'y' in input(f"Overwrite {path} (y/n):"):
                plt.savefig(f'/Users/burggraf1/PycharmProjects/nim2021/nim/manuscript/figs/{path}', dpi=300)

    if plt.rcParams['text.usetex']:
        plt.subplots_adjust(bottom=0.25)


    plt.show()

