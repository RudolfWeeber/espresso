from __future__ import print_function

from math import cos, pi, sin
import numpy as np
import os
import sys

import espressomd
from espressomd import assert_features, lb
from espressomd.lbboundaries import LBBoundary
from espressomd.shapes import Cylinder, Wall, HollowCone


assert_features(["LB_GPU", "LB_BOUNDARIES_GPU"])

# Setup constants

outdir = "./RESULTS_RECTIFICATION_GEOMETRY/"
try:
    os.makedirs(outdir)
except:
    print("INFO: Directory \"{}\" exists".format(outdir))

# Setup the box (we pad the diameter to ensure that the LB boundaries
# and therefore the constraints, are away from the edge of the box)

length = 100
diameter = 20
dt = 0.01

# Setup the MD parameters

system = espressomd.System(box_l=[length, diameter + 4, diameter + 4])
system.cell_system.skin = 0.1
system.time_step = dt
system.min_global_cut = 0.5

# Setup LB parameters (these are irrelevant here) and fluid

agrid = 1
vskin = 0.1
frict = 20.0
visco = 1.0
densi = 1.0

lbf = lb.LBFluidGPU(agrid=agrid, dens=densi, visc=visco, tau=dt, fric=frict)
system.actors.add(lbf)

################################################################################
#
# Now we set up the three LB boundaries that form the rectifying geometry.
# The cylinder boundary/constraint is actually already capped, but we put
# in two planes for safety's sake. If you want to create an cylinder of
# 'infinite length' using the periodic boundaries, then the cylinder must
# extend over the boundary.
#
##########################################################################

# Setup cylinder

cylinder = LBBoundary(
    shape=Cylinder(
        center=[length / 2.0, (diameter + 4) / 2.0, (diameter + 4) / 2.0],
                                     axis=[1, 0, 0],
                                     radius=diameter / 2.0,
                                     length=length,
                                     direction=-1))
system.lbboundaries.add(cylinder)

# Setup walls

wall = LBBoundary(shape=Wall(dist=2, normal=[1, 0, 0]))
system.lbboundaries.add(wall)

wall = LBBoundary(shape=Wall(dist=-(length - 2), normal=[-1, 0, 0]))
system.lbboundaries.add(wall)

# Setup cone

irad = 4.0
angle = pi / 4.0
orad = (diameter - irad) / sin(angle)
shift = 0.25 * orad * cos(angle)

hollow_cone = HollowCone(
    center=(length / 2.0 - shift, (diameter + 4) / 2.0, (diameter + 4) / 2.0),
    axis=[-1, 0, 0],
    outer_radius=orad,
    inner_radius=irad,
    width=2.0,
    opening_angle=angle,
    direction=1)
system.lbboundaries.add(LBBoundary(shape=hollow_cone))

##########################################################################

# Output the geometry

lbf.print_vtk_boundary("{}/boundary.vtk".format(outdir))

##########################################################################
