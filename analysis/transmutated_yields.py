import matplotlib.pyplot as plt
from JSB_tools import decay_nuclide
from JSB_tools import Nuclide, FissionYields
from pathlib import Path
from JSB_tools.MCNP_helper.outp_reader import OutP
import scipy
import numpy as np
from uncertainties import unumpy as unp
"""
Nucei whos rate matches that predicted here:
La-144 - 397 KeV

"""
c_per_second = (192/3.0)*1E-6
charge_per_electron = 1.602E-19
n_electrons = 3*c_per_second / charge_per_electron


ff = 'Sr94'   #  Kr89
decay_rate = True
fraction_escape = 0.1

ff = Nuclide.from_symbol(ff)

outp = OutP(Path(__file__).parent.parent/'mcnp'/'sims'/'du_shot131'/'outp')
tally_down = outp.get_f4_tally('Active down')
u238 = Nuclide.from_symbol('U238')

yields = FissionYields('U238', 'gamma', tally_down.energies)
# yields_cumm = FissionYields('U238', 'gamma', tally_up.energies, independent_bool=False)

# _weights = tally_down.fluxes * u238.gamma_induced_fiss_xs.interp(tally_down.energies)
_weights = tally_down.cell.atom_density * tally_down.dx_per_src * n_electrons * \
           u238.gamma_induced_fiss_xs.interp(tally_down.energies) * fraction_escape

yields.weight_by_erg(_weights)


times = np.arange(0, 300, 1)
data = {}
max_yield = None
i = 0


for parent, yield_ in yields.yields.items():
    yield_ = sum(unp.nominal_values(yield_))

    if max_yield is None:
        max_yield = yield_
    elif yield_ < max_yield*0.01:
        break

    f = decay_nuclide(parent, decay_rate)
    decay_yields = f(times)

    if parent not in data:
        data[parent] = {}

    data[parent]['fission'] = yield_*decay_yields[parent]
    del decay_yields[parent]

    for daughter, v in decay_yields.items():
        if daughter not in data:
            data[daughter] = {}

        try:
            data[daughter][parent] += yield_ * v
        except KeyError:
            data[daughter][parent] = yield_ * v
    i += 1



def srt(dict_):
    try:
        out = {'fission': dict_['fission']}
        del dict_['fission']

    except KeyError:
        out = {}

    for k, v in sorted(dict_.items(), key=lambda k_v: -sum(k_v[1])):
        out[k] = v
    return out


if __name__ == '__main__':
    print(f"Fission fragment: {ff}")
    print(f"Decay lines from FF {ff.name}")
    for g in ff.decay_gamma_lines[:3]:
        print(g)

    fig, ax = plt.subplots()

    y = None

    max_plt_value = max(np.sum([v for v in data[ff.name].values()], axis=0))

    for k, v in srt(data[ff.name]).items():
        if k != 'fission':
            n = Nuclide.from_symbol(k)
            parent_hl = ' ' + n.human_friendly_half_life(include_errors=False)
            p = None
            for g in n.decay_gamma_lines[:3]:
                if g.intensity > 0.01:
                    if p is None:
                        p = True
                        print(n)
                    print(f'\t{g}')
        else:
            parent_hl = ''

        if y is None:
            y = np.zeros_like(v)

        if max(v)/max_plt_value<0.01:
            print(f"Yield from {k} too small. Not included in plot. ")
            continue

        ax.fill_between(times, y, y + v, label=f'{k}-{parent_hl}')
        ax.plot(times, y + v, c='black', linewidth=1)

        y += v

    ax.set_xlabel("Time [s]")

    if decay_rate:
        ax.set_ylabel(f"Decay rate of {ff.name} nuclei [arb. units]")
    else:
        ax.set_ylabel(f"# of {ff.name} nuclei [arb. units]")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(reversed(handles), reversed(labels), title='Source FF')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.show()

