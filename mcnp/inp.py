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
chamber_mount_width = 3/8*inch
dist2chamber_begin = 15/16*inch  # distance from W convertor 2 point where internal cylinder of the chamber begins
chamber_cap_thickness = 0.5
chamber_target_dist = 1.1  # distance of target in chamber from Ti window.
vcd_lead = 8*inch
chamber_target_is_nickel = True
vdc_nickel_radius = 1
beam_fwhm = 0.5
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


if chamber_target_is_nickel:
    chamber_target_radius = 0.6
else:
    chamber_target_radius = 0.5
chamber_target_w = 0.033/(np.pi*chamber_target_radius**2)/nickel_mat.density  # thickness of target in chamber when nickel was the target.

vcd_nickel_dx = 1 * inch
vcd_nickel_dy = 1.5 * inch
vcd_nickel_dz = 0.18/(vcd_nickel_dx*vcd_nickel_dy*nickel_mat.density)

imp = ("ep", 1)

w_convertor = RightCylinder(2, w_mat, imp, z0=0, dz=-0.1*inch, cell_name='Convertor')

source_ccc = RightCylinder(radius=1, z0=w_convertor.z0-2.4*inch, dz=dz,
                           cell_name='CCC source cell', importance=imp)


# ========================= BEGIN CHAMBER GEOM =========================
#  Chamber begins whee Ti begins
chamber = RightCylinder(2, z0=dist2chamber_begin, dz=eff_chamber_length, material=gas_mat, importance=imp,
                        cell_name='Chamber (gas)')

aperture_dummy = RightCylinderSurface.from_min_max_coords(0.5, z0=w_convertor.z0, z1=chamber.z1 + 1,
                                                          surf_name='Aperture dummy')

chamber_mount_up = RightCylinder(6, z0=chamber.z0, dz=chamber_mount_width, material=steel_mat, importance=imp,
                                 cell_name='Mount up')
chamber_mount_up.geometry = -chamber_mount_up & +chamber.surface
chamber_mount_down = RightCylinder(6, z0=chamber.z1, dz=-chamber_mount_width, material=steel_mat, importance=imp,
                                   cell_name='Mount down')
chamber_mount_down.geometry = -chamber_mount_down & +chamber.surface

cap_up = RightCylinder(2, material=steel_mat, importance=imp, z0=chamber.z0, dz=-chamber_cap_thickness,
                       cell_name='Cap up')
cap_down = RightCylinder(2, material=steel_mat, importance=imp, z0=chamber.z1, dz=chamber_cap_thickness,
                         cell_name='Cap down')
cap_up.geometry = -cap_up & +aperture_dummy
cap_down.geometry = -cap_down & +aperture_dummy

ti_up = RightCylinder(0.5, z0=cap_up.z_mid, dz=15*units.um, material=ti_mat, importance=imp, cell_name='Ti up')
ti_down = RightCylinder(0.5, z0=cap_down.z_mid, dz=15*units.um, material=ti_mat, importance=imp,
                        cell_name='Ti down')

chamber_target = RightCylinder(chamber_target_radius, material=nickel_mat, importance=imp, z0=chamber.z0 + chamber_target_dist,
                               dz=chamber_target_w, cell_name='Chamber target')
target_mount = CuboidCell(-1, +1, -1, 1, z0=chamber_target.z0, dz=-0.05, importance=imp, material=steel_mat,
                          cell_name='Target mount')
target_mount.geometry = -target_mount.surface & +aperture_dummy.surface

chamber.geometry = -chamber & +chamber_target & ~target_mount  # set chamber geom now that target is added
# ========================= END CHAMBER GEOM =========================

vcd_lead = RightCylinder(room_w, material=pb_mat, importance=('ep', 20), z0=w_convertor.z1 + 102, dz=vcd_lead,
                         cell_name='VCD lead')

vcd_nickel = CuboidCell(-vcd_nickel_dx/2, vcd_nickel_dx/2, -vcd_nickel_dy/2, vcd_nickel_dy/2,
                        z0=vcd_lead.z0, dz=-vcd_nickel_dz,
                        cell_name='VCD nickel', importance=imp, material=nickel_mat)

vcd_cell = RightCylinder(6, z0=137, dz=1.0, importance=('ep', 20), cell_name='VCD cell')

room = RightCylinder.from_min_max_coords(radius=room_w, z0=w_convertor.z0-2.4*inch, z1=vcd_lead.z0,
                                         material=air_mat, cell_name='Room', importance=imp)
room.geometry = -room.surface & +w_convertor.surface & +chamber & ~chamber_mount_up & ~chamber_mount_down & ~cap_down & ~cap_up & \
                ~ti_down & ~ti_up & +vcd_nickel.surface & +source_ccc.surface

room_2 = RightCylinder.from_min_max_coords(vcd_lead.radius, z0=vcd_lead.z1, z1=Surface.global_zmax(), material=air_mat,
                                           importance=('ep', 20), cell_name='Room 2')
room_2.geometry = -room_2.surface & +vcd_cell.surface

universe = Cell(geometry=+room.surface & +room_2.surface & +vcd_lead.surface,
                importance=("ep", 0), cell_name='Universe')


tally_chamber_target = F4Tally(chamber_target, particle='p', tally_name='Chamber target')
tally_chamber_target.set_erg_bins(1, 22, 20)

tally_vcd_nickel = F4Tally(vcd_nickel, particle='p', tally_name='VCD nickel')
tally_vcd_nickel.set_erg_bins(1, 21.5, 20)

tally_vcd_cell = F4Tally(vcd_cell, particle='p', tally_name='VCD cell')
tally_vcd_cell.set_erg_bins(0.1, 10, 15)

CylFMESH('p', rmaxs=(2, room_w), axis_lengths=(10, Surface.global_zmax()), rbins=(15, 35), axis_bins=(10, 50),
         origin=(0, 0, source_ccc.z0))
CylFMESH('e', rmaxs=(2, room_w), axis_lengths=(10, Surface.global_zmax()), rbins=(15, 35), axis_bins=(10, 50),
         origin=(0, 0, source_ccc.z0))
if __name__ == '__main__':
    inp = InputDeck.mcnp_input_deck(Path.cwd()/'inp', new_file_dir=Path.cwd()/'sims')
    inp.write_inp_in_scope(globals())









