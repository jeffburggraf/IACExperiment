import pickle

import numpy as np
from scipy.integrate import quad
from JSB_tools.nuke_data_tools.gamma_spec import gamma_search
from JSB_tools.nuke_data_tools import FissionYields, Nuclide
from JSB_tools.MCNP_helper import OutP
from pathlib import Path
from uncertainties import unumpy as unp
from JSB_tools import decay_nuclide


# ====================================
#  ====================================

outp = OutP(Path(__file__).parent.parent/'mcnp'/'sims'/'du_shot134'/'outp')
tally = outp.get_f4_tally('Active down')


xs = unp.nominal_values(Nuclide.from_symbol('U238').gamma_induced_fiss_xs.interp(tally.energies))
weights = xs*tally.nominal_fluxes


y = FissionYields('U238', 'gamma', tally.energies)
y.weight_by_erg(weights)
y.threshold()

times = np.arange(10, 400, 1)
rates = {}


predicted_time_dep = {}


for name, yields in y.yields.items():
    yield_ = sum(unp.nominal_values(yields))
    for n, rel_hz in decay_nuclide(name, True)(times).items():
        hz = yield_ * rel_hz
        value = np.trapz(hz, times)

        try:
            rates[n] += value
            predicted_time_dep[n] += hz
        except KeyError:
            rates[n] = value
            predicted_time_dep[n] = hz


for n in list(rates.keys()):
    gammas = Nuclide.from_symbol(n).decay_gamma_lines
    if not len(gammas):
        del rates[n]
        continue
    else:
        rates[n] *= gammas[0].intensity.n


rates = {k: v for k, v in filter(lambda k_v: k_v[1] != 0, rates.items())}
rates = {k: v for k, v in sorted(rates.items(), key=lambda k_v: -k_v[1])}


for k in predicted_time_dep.keys():
    predicted_time_dep[k] /= max(predicted_time_dep[k])


with open(Path(__file__).parent/'predicted_time_dep.pickle', 'wb') as f:
    pickle.dump(predicted_time_dep, f)
    pickle.dump(times, f)


all_gammas = []
all_gammas_yields = []


for k, v in rates.items():
    n = Nuclide.from_symbol(k)
    gammas = n.decay_gamma_lines[:3]
    print(f"{n} - {v:.2e}")
    for g in gammas:
        print(f"\t{g}")
        all_gammas.append(g)
        all_gammas_yields.append(v)


all_gammas_yields = np.array(all_gammas_yields)
all_gammas_yields /= max(all_gammas_yields)
gamma_ergs = [g.erg.n for g in all_gammas]
srt = np.argsort(gamma_ergs)

for i in srt:
    g = all_gammas[i]
    y = all_gammas_yields[i]
    if y < 1E-5:
        continue
    # try:
    hl = Nuclide.from_symbol(g.parent_nuclide_name).half_life.n
    # except AttributeError:
    #     hl = -1
    print(f'Rel. yield: {y:.1e}; t-1/2 = {hl:<8}; {g}')
