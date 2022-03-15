import pickle

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from JSB_tools.nuke_data_tools.gamma_spec import gamma_search
from JSB_tools.nuke_data_tools import FissionYields, Nuclide
from JSB_tools.MCNP_helper import OutP
from pathlib import Path
from uncertainties import unumpy as unp
from JSB_tools import decay_nuclide
from JSB_tools.spectra import EfficiencyCalMixin

_p = Path('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/cal_files/shot120.eff')
effs = EfficiencyCalMixin.stand_alone(np.linspace(0, 3000, 3000), _p)

# effs.plot_efficiency()
# plt.show()

# ====================================
tmin = 0
tmax = 100
ntime_bins = 40
srt_erg = True  # If True sort by gamma energy, else by yield.
gamma_thresh = 5E-3
#  ====================================
print(f"Timerange: {tmin} - {tmax}|")

outp = OutP(Path(__file__).parent.parent/'mcnp'/'sims'/'du_shot134'/'outp')
tally = outp.get_f4_tally('Active down')


xs = unp.nominal_values(Nuclide.from_symbol('U238').gamma_induced_fiss_xs.interp(tally.energies))
weights = xs*tally.nominal_fluxes


y = FissionYields('U238', 'gamma', tally.energies)
y.weight_by_erg(weights)
y.threshold(frac_of_max=0.01)

times = np.linspace(tmin, tmax, ntime_bins)
n_decays = {}


predicted_time_dep = {}


for name, yields in y.yields.items():
    yield_ = sum(unp.nominal_values(yields))
    for n, rel_hz in decay_nuclide(name, 1)(times).items():
        hz = yield_ * rel_hz
        tot_decays = np.trapz(hz, times)

        try:
            n_decays[n] += tot_decays
            predicted_time_dep[n] += hz
        except KeyError:
            n_decays[n] = tot_decays
            predicted_time_dep[n] = hz


# for n in list(rates.keys()):
#     gammas = Nuclide.from_symbol(n).decay_gamma_lines
#     if not len(gammas):
#         del rates[n]
#         continue
#     else:
#         rates[n] *= gammas[0].intensity.n


n_decays = {k: v for k, v in filter(lambda k_v: k_v[1] != 0, n_decays.items())}
n_decays = {k: v for k, v in sorted(n_decays.items(), key=lambda k_v: -k_v[1])}
max_rate = max(n_decays.values())
n_decays = {k: v / max_rate for k, v in n_decays.items()}

for k in predicted_time_dep.keys():
    predicted_time_dep[k] /= max(predicted_time_dep[k])


with open(Path(__file__).parent/'predicted_time_dep.pickle', 'wb') as f:
    pickle.dump(predicted_time_dep, f)
    pickle.dump(times, f)


all_gammas = []
all_gammas_yields = []


for k, v in n_decays.items():
    n = Nuclide.from_symbol(k)
    gammas = filter(lambda x: x.erg.n > gamma_thresh, n.decay_gamma_lines)
    # print(f"{n} - {v:.2e}")

    for g in gammas:
        # print(f"\t{g}")
        all_gammas.append(g)
        all_gammas_yields.append(v*g.intensity*effs(g.erg.n).n)


all_gammas_yields = np.array(all_gammas_yields)
all_gammas_yields /= max(all_gammas_yields)

mask = all_gammas_yields>5E-3
all_gammas = np.array(all_gammas)[mask]
all_gammas_yields = all_gammas_yields[mask]

if srt_erg:
    gamma_ergs = [g.erg.n for g in all_gammas]
    srt = np.argsort(gamma_ergs)
else:
    srt = np.argsort(-all_gammas_yields)


for i in srt:
    g = all_gammas[i]
    y = all_gammas_yields[i]

    number_decay = n_decays[g.parent_nuclide_name]
    if y < 1E-5:
        continue
    # try:
    hl = Nuclide.from_symbol(g.parent_nuclide_name).half_life.n
    # except AttributeError:
    #     hl = -1
    print(f'γ yield: {y:.1e}; decay yield: {number_decay:.1e}; t-1/2 = {hl:<8}; {g}')
