from JSB_tools import Nuclide, FissionYields
from pathlib import Path
from JSB_tools.MCNP_helper.outp_reader import OutP
import scipy

outp = OutP(Path(__file__).parent/'sims'/'1_inp'/'outp')
tally_up = outp.get_tally('Active up')
u238 = Nuclide.from_symbol('U238')

yields = FissionYields('U238', 'gamma', tally_up.energies)
yields_cumm = FissionYields('U238', 'gamma', tally_up.energies, independent_bool=False)

_weights = tally_up.fluxes*u238.gamma_induced_fiss_xs.interp(tally_up.energies)
_weights /= sum(_weights)
yields.weight_by_erg(_weights)
yields_cumm.weight_by_erg(_weights)

for n_name, yield_ in yields.yields.items():
    n = Nuclide.from_symbol(n_name)
    outs = []
    if not 15 <= n.half_life <= 60*4:
        continue
    for g in n.decay_gamma_lines:
        if 218 <= g.erg.n <= 222 and g.intensity.n > 0.01:
            outs.append(g)
    if len(outs):
        print(n, 'Ratio: ', sum(yield_)/sum(yields.yields['Xe139']))
        print(outs)


n = Nuclide.from_symbol('Xe139')
parent_yield = sum(yields[n.name])
print(n)
print(f'cum {n.name} yield: {sum(yields_cumm[n.name])}  indp. {n.name} yield: {sum(yields[n.name])}')

for x in n.get_decay_parents():
    msg = f'{x}  parent yield/daughter yield: {(sum(yields[x.name]))/parent_yield},  '

    print(msg)