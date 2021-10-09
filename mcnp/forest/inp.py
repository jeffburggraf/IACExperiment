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

foil = CuboidCell(-0.5, 0.5, -0.5, 0.5, z0=0, z1=0.05, material=du_mat, importance=imp)

det_cell = CuboidCell(x0=-det_heigth / 2., x1=det_heigth / 2, y0=det_distance, y1=det_distance + det_depth,
                      z0=-det_width / 2, z1=det_width / 2, importance=imp, cell_name='Detector')

room_cell = CuboidCell(x0=det_cell.xmin, x1=det_cell.xmax, y0=foil.ymin, y1=det_cell.ymax,
                       z0=det_cell.zmin, z1=det_cell.zmax, cell_name='Room', importance=imp, material=air_mat)
room_cell.geometry = +foil & +det_cell.surface & -room_cell.surface

universe = Cell(geometry=+room_cell, importance=('np', 0))

tally_n = F4Tally(det_cell, 'n', )
tally_n.set_erg_bins(erg_bins_array=np.arange(0, 9, 0.5))

tally_p = F4Tally(det_cell, 'p')
tally_p.set_erg_bins(erg_bins_array=np.arange(0, 9, 0.5))

if __name__ == '__main__':
    i = InputDeck.mcnp_input_deck(Path.cwd()/'inp')
    i.write_inp_in_scope(globals())