import matplotlib.pyplot as plt

from JSB_tools import Nuclide, FissionYields
from pathlib import Path
from JSB_tools.MCNP_helper.outp_reader import OutP
import scipy

outp = OutP(Path(__file__).parent/'sims'/'du_shot131'/'outp')
tally_up = outp.get_f4_tally('Active down')
u238 = Nuclide.from_symbol('U238')

yields = FissionYields('U238', 'gamma', tally_up.energies)
# yields_cumm = FissionYields('U238', 'gamma', tally_up.energies, independent_bool=False)

_weights = tally_up.fluxes*u238.gamma_induced_fiss_xs.interp(tally_up.energies)
# _weights /= sum(_weights)
yields.weight_by_erg(_weights)
# yields_cumm.weight_by_erg(_weights)

yields_ = {}
for n, yield_ in yields.yields.items():
    yield_ = sum(yield_)
    n = Nuclide.from_symbol(n)
    fac = 0.5**(10/n.half_life)-0.5**(300/n.half_life)

    l = [g.intensity for g in n.decay_gamma_lines if 60 <= g.erg <= 1450]

    if len(l):
        fac *= max(l)
    else:
        continue

    yields_[n.name] = yield_*fac

yields = {k: v for k, v in sorted(yields_.items(), key=lambda x: -x[1])}

i=0

zdata = {}
adata = {}

for k, v in yields.items():
    n = Nuclide.from_symbol(k)
    try:
        zdata[n.Z] += v.n
    except KeyError:
        zdata[n.Z] = v.n
    try:
        adata[n.A] += v.n
    except KeyError:
        adata[n.A] = v.n

    print(f'{k}:{v}; Z={n.Z}, {n}')
    i += 1
    if i > 20:
        break

plt.bar(list(zdata.keys()), list(zdata.values()))
plt.figure()
plt.bar(list(adata.keys()), list(adata.values()))
plt.show()