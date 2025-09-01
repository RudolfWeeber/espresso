#
# Copyright (C) 2013-2022 The ESPResSo project
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
Demonstrates the ESPREsSo to ASE interface
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

box_l = 10.7437
density = 0.7

# Interaction parameters (repulsive Lennard-Jones)

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 2.**(1/6) * lj_sig

# Integration parameters
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

system.time_step = 0.005
system.cell_system.skin = 0.4

warm_steps=10
int_steps = 100
int_n_times = 50




# Particle setup

volume = system.volume()
n_part = int(volume * density)

system.part.add(pos=np.random.random((n_part,3)) * system.box_l)


## Setup of ase interface and Lennard-Jones calculator

# Mapping of ESPResSo types to ASE types

type_mapping = {0:0} 

ase = espressomd.plugins.ase.ASEInterface(
     system, type_mapping, system.part.all())

# ASE calculator tor provide Lennard-Jones forces
lj = LennardJones(sigma=lj_sig,    
                    epsilon=lj_eps,    
                    rc=lj_cut,        
                    smooth=False) 
ase.atoms.calc = lj
# Overlap removal via steepest descent

system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                       max_displacement=lj_sig / 100)



# Integraet one step to have valid particle forces
ase.integrate(1,lj)

max_force =  np.abs(np.amax(system.part.all().f))
while max_force > 10:
    ase.integrate(warm_steps,lj)
    max_force =  np.abs(np.amax(system.part.all().f))
    print("Overlap removal:",max_force)


# Normal integration


system.thermostat.set_langevin(kT=1, gamma=0.2, seed=42)
system.integrator.set_symplectic_euler()

for i in range(int_n_times):
    print(f"run {i} at time={system.time:.2f}")

    ase.integrate(int_steps,lj)
    print(1/n_part *np.sum(system.part.all().v**2)/2)
