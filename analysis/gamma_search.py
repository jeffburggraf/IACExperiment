import numpy as np
from JSB_tools.nuke_data_tools import Nuclide, DecayNuclide
from JSB_tools.nuke_data_tools.nuclide.fission_yields import FissionYields
from JSB_tools.MCNP_helper import OutP
from pathlib import Path
from uncertainties import unumpy as unp
from JSB_tools.spectra import EfficiencyCalMixin

_p = Path('/Users/burggraf1/PycharmProjects/IACExperiment/efficiencies/eff_main.eff')
effs = EfficiencyCalMixin()
effs.unpickle_efficiency(_p)

"""

La144, 397 keV
  4.5 +/- 0.1
  4.8 +/- 0.2
  3.5 +/- 0.2
  3.2 +/- 0.2

Good consistency for Sr93 @ 590 keV


"""

# effs.plot_efficiency()
# plt.show()

c_per_second = (192/3.0)*1E-6
charge_per_electron = 1.602E-19
n_electrons = 3*c_per_second / charge_per_electron

# ====================================
tmin = 100
tmax = 300
ntime_bins = 200
srt_erg = True  # If True sort by gamma energy, else by yield.
gamma_erg_range = 50, 2000
rel_yield_thresh = 3E-3
eff_corr = True
#  ====================================
print(f"Time range: {tmin} - {tmax}|")


outp = OutP(Path(__file__).parent.parent/'mcnp'/'sims'/'du_shot134'/'outp')
tally = outp.get_f4_tally('Active down')


xs = unp.nominal_values(Nuclide('U238').gamma_induced_fiss_xs(tally.energies))
weights = xs*tally.dx_per_src*n_electrons*tally.cell.atom_density*0.1


y = FissionYields('U238', 'gamma', tally.energies)
y.weight_by_erg(weights)
# y.threshold(frac_of_max=0.001)

times = np.linspace(tmin, tmax, ntime_bins)
n_decays = {}

for n, yield_ in y.yields.items():
    yield_ = sum(unp.nominal_values(yield_))
    for n, hz in DecayNuclide(n)(times, decay_rate=True).items():
        tot_decays = yield_*np.trapz(hz, times)

        try:
            n_decays[n] += tot_decays
        except KeyError:
            n_decays[n] = tot_decays

gamma_lines = []
gamma_yields = []

for n, decays in n_decays.items():
    nuclide = Nuclide(n)
    gs = nuclide.decay_gamma_lines
    for g in nuclide.decay_gamma_lines:
        if not gamma_erg_range[0] <= g.erg.n <= gamma_erg_range[1]:
            continue
        y = g.intensity*decays
        if eff_corr:
            y *= effs(g.erg.n)
        gamma_yields.append(y)
        gamma_lines.append(g)

if not srt_erg:
    srt = np.argsort(gamma_yields)[::-1]
else:
    srt = np.argsort([g.erg.n for g in gamma_lines])

max_yield = max(gamma_yields)


for i in srt:
    g = gamma_lines[i]
    abs_g_yield = gamma_yields[i]
    rel_g_yield = abs_g_yield/max_yield

    if rel_g_yield < rel_yield_thresh:
        continue

    number_decay = n_decays[g.parent_nuclide_name]

    hl = Nuclide(g.parent_nuclide_name).half_life.n
    erg = f"{g.erg:.2f} keV"
    print(f'{erg: <9}; rel/abs γ yield: {rel_g_yield:.1e} / {abs_g_yield.n:.1e}; t-1/2 = {hl:<8}; '
          f'decay yield: {number_decay:.1e};  {g}')

