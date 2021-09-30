import numpy as np
from JSB_tools.MCNP_helper.materials import DepletedUranium, Material, Titanium, Tungsten, StainlessSteel, Air, Mylar
from JSB_tools.MCNP_helper.geometry.primitives import RightCylinder, RightCylinderSurface, SphereSurface
from JSB_tools.MCNP_helper import F4Tally, units, CuboidCell, Cell, Surface
from JSB_tools.MCNP_helper.input_deck import MCNPSICard, InputDeck
from pathlib import Path
from JSB_tools.MCNP_helper.units import inch
from erg_dist import neutron_ergs, neutron_fluxes

#  ====================================
det_distance = 125
#  ====================================
det_width = 6*inch
det_heigth = 30*inch
det_depth = 1.5*inch

si_card = MCNPSICard(neutron_ergs, neutron_fluxes/np.sum(neutron_fluxes))

imp = ('np', 1)

du_mat = DepletedUranium()
air_mat = Air()

foil = CuboidCell(-0.5, 0.5, -0.5, 0.5, zmin=0, zmax=0.05, material=du_mat, importance=imp)

det_cell = CuboidCell(xmin=-det_heigth/2., xmax=det_heigth/2, ymin=det_distance, ymax=det_distance+det_depth,
                      zmin=-det_width/2, zmax=det_width/2, importance=imp, cell_name='Detector')

room_cell = CuboidCell(xmin=det_cell.xmin, xmax=det_cell.xmax, ymin=foil.ymin, ymax=det_cell.ymax,
                       zmin=det_cell.zmin, zmax=det_cell.zmax, cell_name='Room', importance=imp, material=air_mat)
room_cell.geometry = +foil & +det_cell.surface & -room_cell.surface

universe = Cell(geometry=+room_cell, importance=('np', 0))

tally_n = F4Tally(det_cell, 'n', )
tally_n.set_erg_bins(erg_bins_array=np.arange(0, 9, 0.5))

tally_p = F4Tally(det_cell, 'p')
tally_p.set_erg_bins(erg_bins_array=np.arange(0, 9, 0.5))

if __name__ == '__main__':
    i = InputDeck.mcnp_input_deck(Path.cwd()/'inp')
    i.write_inp_in_scope(globals())