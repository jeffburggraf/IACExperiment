import numpy as np
from JSB_tools.MCNP_helper.input_deck import InputDeck
from JSB_tools.MCNP_helper.geometry.geom_core import Cell, Surface
from JSB_tools.MCNP_helper.materials import (DepletedUranium, Material, Titanium, Tungsten, StainlessSteel, Air,
                                             Lead, Nickel, Aluminum)
from JSB_tools.MCNP_helper.geometry.primitives import RightCylinder, RightCylinderSurface
from JSB_tools.MCNP_helper import F4Tally, units, CellGroup, CuboidCell, CylFMESH, TRCL, CuboidSurface
from pathlib import Path
from JSB_tools.MCNP_helper.units import inch

# ==========================================
w_to_target_dist = 1*inch + 1.1
src_to_w_distance = 2.4*inch
# ==========================================
dz=1E-4
imp = ('ep', 1)


w_mat = Tungsten()
air_mat = Air()
ni_mat = Nickel()
du_mat = DepletedUranium()
ti_mat = Titanium()

ni_target_w = 0.005*inch
ni_target_radius = np.sqrt(0.033/(np.pi*ni_target_w*ni_mat.density))

convertor = RightCylinder(2, z0=src_to_w_distance, dz=0.1*inch, material=w_mat, importance=imp,
                          cell_name='Convertor')
ti_window = RightCylinder(2, z0=convertor.z1 + 0.75*inch, dz=15*units.um, material=ti_mat, importance=imp,
                          cell_name='Ti')
room = RightCylinder(2, z0=-dz, material=air_mat, importance=imp, cell_name='Room')
room.dz = w_to_target_dist + src_to_w_distance + 1

target = RightCylinder(0.5, z0=convertor.z1 + w_to_target_dist, dz=0.05,  material=du_mat, importance=imp,
                       cell_name='Target')
room.geometry = -room & + target & +convertor & +ti_window

world = Cell(geometry=+room, importance=("pe", 0))

tally = F4Tally(target, 'p', tally_name='Target')
tally.set_erg_bins(4, 22, 20)

input_deck = InputDeck.mcnp_input_deck(Path(__file__).parent/'inp_simple',
                                       new_file_dir=Path(__file__).parent/'sims')
input_deck.write_inp_in_scope(globals(), 'simple_du')

target.radius = ni_target_radius
target.dz = ni_target_w
target.material = ni_mat
input_deck.write_inp_in_scope(globals(), 'simple_ni')

