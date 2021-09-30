"""
MCNP model of chamber and related during the IAC experiment
"""
import numpy as np
from JSB_tools.MCNP_helper.input_deck import InputDeck
from JSB_tools.MCNP_helper.geometry.geom_core import Cell, Surface
from JSB_tools.MCNP_helper.materials import DepletedUranium, Material, Titanium, Tungsten, StainlessSteel, Air, Lead, Nickel
from JSB_tools.MCNP_helper.geometry.primitives import RightCylinder, RightCylinderSurface
from JSB_tools.MCNP_helper import F4Tally, units, CellGroup, CuboidCell, CylFMESH, TRCL
from pathlib import Path
from JSB_tools.MCNP_helper.units import inch

# ======================================================
dz=1E-4
room_w = 6*inch
o_ring_width = 1/8*inch  # Width of o-rings
chamber_mount_width = 1
W2Ti_distance = 3/4*inch  # distance from W convertor to Ti window.
chamber_cap_thickness = 0.5
chamber_target_dist = 1.1  # distance of target in chamber from Ti window.
vcd_lead = 8*inch
# ======================================================
#  --- Calculated parameters ---
eff_chamber_length = 6 + 4*o_ring_width + 2*chamber_mount_width  # effective chamber length
electron_erg = 21.5
# /---

#  --- MATERIALS ---
w_mat = Tungsten()
# du_mat = DepletedUranium()
steel_mat = StainlessSteel()
air_mat = Air()
ti_mat = Titanium()
gas_mat = Material.gas(['He', 'Ar'], atom_fractions=[0.5, 0.5], pressure=1.16, mat_name='He_Ar')
nickel_mat = Nickel()
pb_mat = Lead()
# /---


chamber_target_w = 0.033/(np.pi*0.5**2)/nickel_mat.density  # thickness of target in chamber when nickel was the target.
vcd_nickel_w = 0.18/(np.pi*1**2)/nickel_mat.density  # radius is 1 cm (this is approx)


imp = ("ep", 1)

w_convertor = RightCylinder(2, w_mat, imp, z0=0, dz=-0.1*inch, cell_name='Convertor')

#  Chamber begins whee Ti begins
chamber = RightCylinder(2, z0=W2Ti_distance, dz=eff_chamber_length, material=gas_mat, importance=imp,
                        cell_name='Chamber (gas)')

aperture_dummy = RightCylinderSurface.from_min_max_coords(0.5, z0=w_convertor.zmin, z1=chamber.zmax + 1,
                                                          surf_name='Aperture dummy')

mount_up = RightCylinder(6, z0=chamber.zmin, dz=chamber_mount_width, material=steel_mat, importance=imp,
                         cell_name='Mount up')
mount_up.geometry = -mount_up & +chamber.surface
mount_down = RightCylinder(6, z0=chamber.zmax, dz=-chamber_mount_width, material=steel_mat, importance=imp,
                           cell_name='Mount down')
mount_down.geometry = -mount_down & +chamber.surface

cap_down = RightCylinder(2, material=steel_mat, importance=imp, z0=chamber.zmax, dz=chamber_cap_thickness,
                         cell_name='Cap down')
cap_down.geometry = -cap_down & +aperture_dummy

cap_up = RightCylinder(2, material=steel_mat, importance=imp, z0=chamber.zmin, dz=-chamber_cap_thickness,
                       cell_name='Cap up')
cap_up.geometry = -cap_up & +aperture_dummy

ti_up = RightCylinder(0.5, z0=chamber.zmin, dz=-15*units.um, material=ti_mat, importance=imp, cell_name='Ti up')
ti_down = RightCylinder(0.5, z0=chamber.zmax, dz=15*units.um, material=ti_mat, importance=('ep', 10),
                        cell_name='Ti down')

chamber_target = RightCylinder(0.5, material=nickel_mat, importance=imp, z0=chamber.zmin + chamber_target_dist,
                               dz=chamber_target_w, cell_name='Chamber target')

chamber.geometry = -chamber & +chamber_target  # set chamber geom now that target is added

vcd_lead = RightCylinder(room_w-1E-3, material=pb_mat, importance=('ep', 20), z0=w_convertor.zmax + 102, dz=vcd_lead,
                         cell_name='VCD lead')

vcd_nickel = RightCylinder(1, material=nickel_mat, importance=imp, z0=vcd_lead.zmin - 1E-2, dz=-vcd_nickel_w,
                           cell_name='VCD nickel')

vcd_cell = RightCylinder(6, z0=137, dz=1.0, importance=imp, cell_name='VCD cell')

room = CuboidCell(-room_w, room_w, -room_w, room_w, zmin=w_convertor.zmin-2*dz, zmax=Surface.global_zmax(),
                  material=air_mat, cell_name='Room', importance=imp)
room.geometry = -room.surface & +w_convertor.surface & +chamber & ~mount_up & ~mount_down & ~cap_down & ~cap_up &\
                ~ti_down & ~ti_up & +vcd_lead.surface & +vcd_nickel.surface & ~vcd_cell.surface

universe = Cell(geometry=+room.surface, importance=("ep", 0), cell_name='Universe')


tally_chamber_target = F4Tally(chamber_target, particle='p', tally_name='Chamber target')
tally_chamber_target.set_erg_bins(1, 21.5, 20)


tally_vcd_nickel = F4Tally(vcd_nickel, particle='p', tally_name='VCD nickel')
tally_vcd_nickel.set_erg_bins(1, 21.5, 20)

tally_vcd_cell = F4Tally(vcd_cell, particle='p', tally_name='VCD cell')
tally_vcd_cell.set_erg_bins(0.1, 10, 15)

CylFMESH('p', 10, Surface.global_zmax(), rbins=20, axis_bins=50)
CylFMESH('e', 10, Surface.global_zmax(), rbins=20, axis_bins=50)

if __name__ == '__main__':
    inp = InputDeck.mcnp_input_deck(Path.cwd()/'inp', new_file_dir=Path.cwd()/'sims')
    inp.write_inp_in_scope(globals())









