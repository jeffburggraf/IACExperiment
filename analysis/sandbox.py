import re

import matplotlib.pyplot as plt
from JSB_tools import DecayNuclide
from JSB_tools import Nuclide, FissionYields
from pathlib import Path
from JSB_tools.MCNP_helper.outp_reader import OutP
import scipy
import numpy as np
import pickle
from uncertainties import unumpy as unp
from JSB_tools import TabPlot



outp = OutP(Path(__file__).parent.parent/'mcnp'/'sims'/'du_shot131'/'outp')
tally_down = outp.get_f4_tally('Active down')
u238 = Nuclide.from_symbol('U238')

yields = FissionYields('U238', 'gamma', tally_down.energies)
# yields_cumm = FissionYields('U238', 'gamma', tally_up.energies, independent_bool=False)
print(yields.library)
# _weights = tally_down.fluxes * u238.gamma_induced_fiss_xs.interp(tally_down.energies)
_weights = tally_down.dx_per_src * \
           u238.gamma_induced_fiss_xs.interp(tally_down.energies)

_weights = unp.nominal_values(_weights/sum(_weights))

yields.weight_by_erg(_weights)


times = np.logspace(0, 5, 100)

max_yield = None

yields.threshold(0.001)

final_rates = {}
# yields_ =

decay_rates = np.zeros_like(times, dtype=float)
n_neutrons = np.zeros_like(decay_rates)

max_tabs = 1
n_tabs = 1
t = TabPlot()

for n, y in yields.yields.items():
    scale = sum(unp.nominal_values(y))
    d = DecayNuclide(n, scale)
    tot_rates = np.zeros_like(decay_rates)
    tot_amounts = np.zeros_like(tot_rates)

    for (k, v), (_, v_decay) in zip(d(times, decay_rate=False).items(), d(times, decay_rate=True).items()):
        ff = Nuclide.from_symbol(k)
        decay_rates += unp.nominal_values(v_decay)
        tot_rates += v_decay
        tot_amounts += v
        N = ff.N
        n_neutrons += N*v

        final_rates[(n, k)] = v_decay[-1]

        if n_tabs <= max_tabs:
            try:
                ax = t.new_ax(ff.latex_name)
            except OverflowError:
                t = TabPlot()
                n_tabs += 1
                ax = t.new_ax(ff.latex_name)

            ax.plot(times, N*v)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_title(ff.latex_name)

    # fig, axs = plt.subplots(1,2)
    # plt.plot(times, tot_amounts, label='Amount')
    # axs[1].plot(times, tot_rates, label='Decay rate')
    # axs[1].legend()
    # axs[1].set_xscale('log')
    # axs[1].set_yscale('log')
    # axs[0].set_xscale('log')
    # axs[0].set_yscale('log')
    # axs[0].plot(times, decay_rates)
    # plt.title(f"{Nuclide.from_symbol(n)} ")
    #
    # plt.show()



final_rates = {k:v for k,v in sorted(final_rates.items(), key=lambda x:-x[-1])}
for k , v in final_rates.items():
    print(k, v)

plt.plot(times, decay_rates)
plt.figure()
plt.plot(times, n_neutrons)

plt.show()