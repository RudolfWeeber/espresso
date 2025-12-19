#
# Copyright (C) 2025 The ESPResSo project
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
#
"""
Demonstrate the ESPResSo to ASE interface.
Simulate a Lennard-Jones fluid maintained at a fixed temperature
by a Langevin thermostat.
"""
import numpy as np
import espressomd
import espressomd.plugins.ase
from ase.calculators.lj import LennardJones

required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

# System parameters
box_l = 10.
lj_eps = 1.0
lj_sig = 1.0
lj_cut = 2.5 * lj_sig

system = espressomd.System(box_l=[box_l] * 3)
system.time_step = 0.005
system.cell_system.skin = 0.4

p1 = system.part.add(pos=[0., 0., 0.])
p2 = system.part.add(pos=[-0.92, 0., 0.])

## Setup of ASE interface and Lennard-Jones calculator

# Mapping of ESPResSo types to ASE types
type_mapping = {0: 0}
ase = espressomd.plugins.ase.ASEInterface(
    system, type_mapping, system.part.all())

# ASE calculator to provide Lennard-Jones forces
lj = LennardJones(sigma=lj_sig, epsilon=lj_eps, rc=lj_cut, smooth=False)
ase.atoms.calc = lj

system.integrator.set_vv()

pos_saved = system.part.all().pos
v_saved = system.part.all().v
ase_forces = []

for i in range(4):
    ase.integrate(1, lj)
    if i != 0:
        ase_forces.append(p1.f)

system.part.all().pos = pos_saved
system.part.all().v = v_saved
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")
es_forces = []
system.part.all().ext_force = [0., 0., 0.]
for i in range(3):
    p1.pos = p1.pos  # force cell structure reconstruction
    system.integrator.run(1, reuse_forces=False)
    es_forces.append(p1.f)

np.testing.assert_allclose(ase_forces, es_forces)
