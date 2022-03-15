import pickle
from JSB_tools.MCNP_helper.input_deck import InputDeck
from JSB_tools.MCNP_helper import Cell, Surface, RightCylinder
from JSB_tools.MCNP_helper.materials import Material, Mylar, DepletedUranium, PHITSOuterVoid
from JSB_tools.PHITS_tools import NucleusSource, CylindricalSource, Nuclide, GaussianEnergyDistribution
from pathlib import Path
from JSB_tools.MCNP_helper.units import um
from JSB_tools import FileManager
from typing import List

#  =================================================================
nuclides = ['Xe139', 'Sr94']
nps = 60000
# dedx_scale = 0.15  # have done 0, 0.15, 0.3, 0.5
dedx_scales = [0.1]
vary_materials = ['gas']  # must contain only the following, ['gas', 'my', 'du']
# which_mat = 'my'
# sim_folder = f'vary{which_mat}{dedx_scale:.1f}'
#  =================================================================

fman = FileManager(recreate=False)

if dedx_scales is not None:
    sim_folder = 'vary' + ''.join([f"{mat}{scale}" for mat, scale in zip(vary_materials, dedx_scales)])
else:
    sim_folder = None

fman.add_path(sim_folder, **dict(zip(vary_materials, dedx_scales)), missing_ok=True)


if dedx_scales is None:
    dedx_scales = vary_materials = []

for mat in ['my', 'gas', 'du']:
    if mat not in vary_materials:
        vary_materials.append(mat)
        dedx_scales.append(1)


du = DepletedUranium()
he_ar = Material.gas(['Ar', 'He'], atom_fractions=[1, 1], pressure=1.1)
mylar_mat = Mylar()

mat_dict = {'my': mylar_mat, 'gas': he_ar, 'du': du}


for mat, dedx_scale in zip(vary_materials, dedx_scales):
    mat = mat_dict[mat]
    mat.set_srim_dedx(scaling=dedx_scale)


u_cell = RightCylinder(0.5, du, z0=1, dz=10*um, cell_name="Du foil")


mylar_cell = RightCylinder(0.5, z0=u_cell.z1, material=mylar_mat, cell_name="Mylar")

chamber = RightCylinder(2, he_ar, z0=0, dz=6, cell_name='Chamber')
chamber.geometry = -chamber & +mylar_cell & +u_cell


void = Cell(material=PHITSOuterVoid(), geometry=+chamber)

with open(Path(__file__).parent/'FF_energies.pickle', 'rb') as f:
    ff_energies = pickle.load(f)


def get_ff_erg(n):
    mean, std = ff_energies[n]['mean'], ff_energies[n]['std']
    return GaussianEnergyDistribution(mean, std)


spacial_dist = CylindricalSource(0.5, z0=u_cell.z0, dz=u_cell.dz)

sources = []
for nuclide in nuclides:
    nuclide = Nuclide.from_symbol(nuclide)
    erg_dist = get_ff_erg(nuclide.name)
    source = NucleusSource(nuclide, erg_dist=erg_dist, spacial_dist=spacial_dist)
    sources.append(source)


sim_dir = Path(__file__).parent
if sim_folder is not None:
    sim_dir /= sim_folder

i = InputDeck.phits_input_deck(Path(__file__).parent/'inp', new_file_dir=sim_dir)

for my_th in [0, 2.5, 5, 10, 20]:
    my_th = my_th*um
    if my_th == 0:
        mylar_cell.disable()
        chamber.geometry = -chamber & +u_cell
    else:
        mylar_cell.enable()
        chamber.geometry = -chamber & +mylar_cell & +u_cell

    mylar_cell.dz = my_th
    new_file_name = i.write_inp_in_scope(globals())
    # fman.add_path(new_file_name.parent/'PTRAC.txt', my_th=my_th)



