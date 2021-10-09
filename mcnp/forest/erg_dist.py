from pathlib import Path
import numpy as np

neutron_ergs = []
neutron_fluxes = []
with open(Path.cwd()/'erg_dist') as f:
    lines = f.readlines()
    for e, v in map(str.split, lines):
        flux = eval(v.replace('^', '**'))
        erg = float(e)
        neutron_fluxes.append(flux)
        neutron_ergs.append(erg)

neutron_fluxes = np.array(neutron_fluxes)