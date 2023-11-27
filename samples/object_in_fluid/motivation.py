# Copyright (C) 2010-2022 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Simulate the motion of flexible red blood cells in a lattice-Boltzmann fluid
with solid obstacles. For more details, see :ref:`Object-in-fluid`.
"""

import espressomd
import espressomd.lb
import espressomd.shapes

required_features = ["WALBERLA", "EXTERNAL_FORCES", "SOFT_SPHERE", "MASS"]
espressomd.assert_features(required_features)

import os
import tqdm
import argparse
import warnings

import object_in_fluid as oif
from object_in_fluid.oif_utils import output_vtk_rhomboid, output_vtk_cylinder

parser = argparse.ArgumentParser()
parser.add_argument("sim", metavar="N", type=int, help="simulation identifier")
args = parser.parse_args()

output_path = os.path.join("output", f"sim{args.sim}")
if not os.path.isdir(output_path):
    os.makedirs(output_path)
    print(f'Saving data to {output_path}')
else:
    warnings.warn(
        f"Folder {output_path} already exists, files will be overwritten")

boxX = 22.0
boxY = 14.0
boxZ = 6.0
time_step = 0.05

system = espressomd.System(box_l=(boxX, boxY, boxZ))
system.time_step = time_step
system.cell_system.skin = 0.2

# creating the template for RBCs
cell_type = oif.OifCellType(
    nodes_file=os.path.join("input", "rbc374nodes.dat"),
    triangles_file=os.path.join("input", "rbc374triangles.dat"),
    system=system, ks=0.04, kb=0.016, kal=0.02, kag=0.9, kv=1.0,
    check_orientation=False, resize=(2.0, 2.0, 2.0))

# creating the RBCs
cell0 = oif.OifCell(cell_type=cell_type,
                    particle_type=0, origin=[5.0, 5.0, 3.0])

# cell-wall interactions
system.non_bonded_inter[cell0.particle_type, 10].soft_sphere.set_params(
    a=0.0001, n=1.2, cutoff=0.1, offset=0.0)

# fluid
lbf = espressomd.lb.LBFluidWalberla(
    agrid=1., density=1., kinematic_viscosity=1.5, tau=system.time_step,
    ext_force_density=[0.025, 0., 0.], single_precision=True)
system.lb = lbf
system.thermostat.set_lb(LB_fluid=lbf, gamma=1.5)

# creating boundaries and obstacles in the channel
# OutputVtk writes a file
# boundaries for the fluid are set up by marking LB nodes as boundaries, here with the help of shapes
# boundaries for the cells are created by creating constraints from the shapes

boundary_shapes = []

# bottom of the channel
bottom_shape = espressomd.shapes.Rhomboid(corner=[0.0, 0.0, 0.0], a=[boxX, 0.0, 0.0],
                                          b=[0.0, boxY, 0.0], c=[0.0, 0.0, 1.0],
                                          direction=1)
boundary_shapes.append(bottom_shape)
output_vtk_rhomboid(
    bottom_shape, out_file=os.path.join(output_path, "wallBottom.vtk"))

# top of the channel
top_shape = espressomd.shapes.Rhomboid(corner=[0.0, 0.0, boxZ - 1], a=[boxX, 0.0, 0.0],
                                       b=[0.0, boxY, 0.0], c=[0.0, 0.0, 1.0], direction=1)
boundary_shapes.append(top_shape)
output_vtk_rhomboid(
    top_shape, out_file=os.path.join(output_path, "wallTop.vtk"))

# front wall of the channel
front_shape = espressomd.shapes.Rhomboid(corner=[0.0, 0.0, 0.0], a=[boxX, 0.0, 0.0],
                                         b=[0.0, 1.0, 0.0], c=[0.0, 0.0, boxZ], direction=1)
boundary_shapes.append(front_shape)
output_vtk_rhomboid(
    front_shape, out_file=os.path.join(output_path, "wallFront.vtk"))

# back wall of the channel
back_shape = espressomd.shapes.Rhomboid(corner=[0.0, boxY - 1.0, 0.0], a=[boxX, 0.0, 0.0],
                                        b=[0.0, 1.0, 0.0], c=[0.0, 0.0, boxZ], direction=1)
boundary_shapes.append(back_shape)
output_vtk_rhomboid(
    back_shape, out_file=os.path.join(output_path, "wallBack.vtk"))

# obstacle - cylinder A
cylA_shape = espressomd.shapes.Cylinder(center=[11.0, 2.0, boxZ / 2.], axis=[0.0, 0.0, 1.0],
                                        length=boxZ, radius=2.0, direction=1)
boundary_shapes.append(cylA_shape)
output_vtk_cylinder(
    cylA_shape, n=20, out_file=os.path.join(output_path, "cylinderA.vtk"))

# obstacle - cylinder B
cylB_shape = espressomd.shapes.Cylinder(center=[16.0, 8.0, boxZ / 2.], axis=[0.0, 0.0, 1.0],
                                        length=boxZ, radius=2.0, direction=1)
boundary_shapes.append(cylB_shape)
output_vtk_cylinder(
    cylB_shape, n=20, out_file=os.path.join(output_path, "cylinderB.vtk"))

# obstacle - cylinder C
cylC_shape = espressomd.shapes.Cylinder(center=[11.0, 12.0, boxZ / 2.], axis=[0.0, 0.0, 1.0],
                                        length=boxZ, radius=2.0, direction=1)
boundary_shapes.append(cylC_shape)
output_vtk_cylinder(
    cylC_shape, n=20, out_file=os.path.join(output_path, "cylinderC.vtk"))

for shape in boundary_shapes:
    lbf.add_boundary_from_shape(shape)
    system.constraints.add(shape=shape, particle_type=10)


def write_cells_vtk(i):
    filepath = os.path.join(output_path, "cell{cell_id}_{index}.vtk")
    cell0.output_vtk_pos_folded(file_name=filepath.format(cell_id=0, index=i))


maxCycle = 100
# main integration loop
for i in tqdm.tqdm(range(maxCycle)):
    write_cells_vtk(i)
    system.integrator.run(steps=100)
write_cells_vtk(maxCycle)
print("Simulation completed.")
