import pickle
from JSB_tools.MCNP_helper.input_deck import InputDeck
from JSB_tools.MCNP_helper import Cell, Surface, RightCylinder
from JSB_tools.MCNP_helper.materials import Material, Mylar, DepletedUranium, PHITSOuterVoid
from JSB_tools.PHITS_tools import NucleusSource, CylindricalSource, Nuclide, GaussianEnergyDistribution
from pathlib import Path
from JSB_tools.MCNP_helper.units import um

#  =================================================================
nuclides = ['Xe139', 'Mo105', ]
nps = 10000
#  =================================================================

du = DepletedUranium()
he_ar = Material.gas(['Ar', 'He'], atom_fractions=[1, 1], pressure=1.1)
mylar_mat = Mylar()
he_ar.set_srim_dedx()
du.set_srim_dedx()
mylar_mat.set_srim_dedx()

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

i = InputDeck.phits_input_deck(Path(__file__).parent/'inp')

for my_th in [0, 2.5, 5, 10, 20]:
    my_th = my_th*um
    if my_th == 0:
        mylar_cell.disable()
        chamber.geometry = -chamber & +u_cell
    else:
        mylar_cell.enable()
        chamber.geometry = -chamber & +mylar_cell & +u_cell

    mylar_cell.dz = my_th
    i.write_inp_in_scope(globals())



