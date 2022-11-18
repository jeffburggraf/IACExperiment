import re
import matplotlib.pyplot as plt
from JSB_tools import TabPlot
from JSB_tools.nuke_data_tools import Nuclide
from analysis import Shot
from uncertainties import unumpy as unp
from pathlib import Path
from JSB_tools.spe_reader import SPEFile


# spe_path = "/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/EffCalFriday/Eu152-0.Spe"
spe_path = "/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/Nickel/Nickel.Spe"

spe_path = Path(spe_path)


if spe_path.name == 'Nickel.Spe':
    n = Nuclide("Ni57")
else:
    n = re.match("([A-Z][a-z]?[0-9]+)-", spe_path.name)
    assert n
    n = Nuclide(n.groups()[0])

eff_path_old = "/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/cal_files/shot121.eff"
eff_path_new = "/Users/burggraf1/PycharmProjects/IACExperiment/efficiencies/angel.eff"


spe = SPEFile(spe_path)
spe.plot_efficiency()


def get_window(erg):
    w0 = 2
    w1 = 3
    return w0 + w1*erg/2000


def omg():
    if line.parent_nuclide.name == 'Co57' and line.erg < 30:
        return False
    return True


xs = []
ys = []

tab = TabPlot()
tab.fig.suptitle(spe_path.name)

axss = []

for line in sorted(n.decay_gamma_lines, key=lambda x: x.erg.n):
    if not omg():
        continue

    if line.intensity > 0.05:
        erg = line.erg.n
        ax1, ax2 = tab.new_ax(f"{erg:.1f}", 1, 2)
        axss.append((ax1, ax2))

        emin = erg - get_window(erg)
        emax = erg + get_window(erg)

        counts = sum(spe.get_counts(erg_min=emin, erg_max=emax, eff_corr=True, remove_baseline=True,
                                    nominal_values=False, debug_plot=ax1, baseline_method='median'))
        ax2.axvline(erg, linewidth=1, ls='--')

        xs.append(erg)
        ys.append(counts/line.intensity.n)


for ax1, ax2 in axss:
    ax2.errorbar(xs, unp.nominal_values(ys), unp.std_devs(ys), ls='None', marker='o')
    ax2.set_ylim(0, ax2.get_ylim()[1])

plt.show()
