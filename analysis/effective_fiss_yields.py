import numpy as np
from matplotlib import pyplot as plt
from JSB_tools import mpl_hist
from pathlib import Path
import pickle
from JSB_tools.nuke_data_tools import FissionYields, Nuclide
from JSB_tools.MCNP_helper import OutP
from pathlib import Path
from uncertainties import unumpy as unp
from uncertainties import ufloat
from JSB_tools import decay_nuclide


thresh_frac = 0.1
outp = OutP(Path(__file__).parent.parent/'mcnp'/'sims'/'du_shot134'/'outp')
tally = outp.get_f4_tally('Active down')

xs = unp.nominal_values(Nuclide.from_symbol('U238').gamma_induced_fiss_xs.interp(tally.energies))
weights = xs*tally.nominal_fluxes
weights /= sum(weights)


y = FissionYields('U238', 'gamma', tally.energies)

y.weight_by_erg(weights)
print("Excluding", y.threshold(0.02))

rates_dict = {}
remaining_dict = {}

dt = 2.5
times = np.arange(20, 300 + dt, dt)

i = 1
entries = len(y.yields)

yields = {k: np.sum(v) for k, v in y.yields.items()}

for n_name, yield_ in yields.items():
    print(f"Calculating {n_name}... {i} out of {entries}")
    i += 1
    nuclide = Nuclide.from_symbol(n_name)
    if nuclide.is_stable or not (1 < nuclide.half_life.n < 15*60):
        continue
    func = decay_nuclide(n_name)

    rates = func(times, yield_, decay_rate=True)

    n_nuclides = func(times[-1], yield_)
    for fp, rate in rates.items():
        try:
            rates_dict[n_name] += rate
            remaining_dict[n_name] += n_nuclides[fp]
        except KeyError:
            rates_dict[n_name] = rate
            remaining_dict[n_name] = n_nuclides[fp]


result = {}
yield_vs_a = {}

i = 1
entries = len(rates_dict)

for k, v in rates_dict.items():
    print(f"Summing {k}... {i} out of {entries}")
    i += 1
    n = Nuclide.from_symbol(k)
    result[k] = np.trapz(v, dx=dt)
    result[k] += remaining_dict[k]
    try:
        yield_vs_a[n.A] += result[k]
    except KeyError:
        yield_vs_a[n.A] = result[k]

max_yield = max(yield_vs_a.values())


yield_vs_a = {k: v for k,v in yield_vs_a.items() if (v >= thresh_frac*max_yield and v.std_dev/v.n<0.3)}



As = list(yield_vs_a.keys())
mass_yield = list(yield_vs_a.values())
mass_yield_err = unp.std_devs(mass_yield)
mass_yield = unp.nominal_values(mass_yield)
plt.errorbar(As, mass_yield, mass_yield_err, ls='None', marker='o')
plt.show()

