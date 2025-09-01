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
import numpy as np
from ase.calculators.calculator import Calculator, all_changes

class HarmonicAllPairs(Calculator):
    """
    Minimal ASE calculator: harmonic springs between all pairs (i < j).

    E = sum_{i<j} 0.5 * k * (r_ij - r0)**2
    F_i = - dE/d r_i = -k * (r_ij - r0) * (r_vec / r_ij)   (equal/opposite on j)

    Parameters
    ----------
    k : float
        Spring constant (energy per distance^2).
    r0 : float
        Rest length for every pair (same for all pairs).
    cutoff : float or None
        If given, ignore pairs with r_ij > cutoff.
    """
    implemented_properties = ['energy', 'forces']

    def __init__(self, k=1.0, r0=1.0, cutoff=None, **kwargs):
        super().__init__(**kwargs)
        self.k = float(k)
        self.r0 = float(r0)
        self.cutoff = float(cutoff) if cutoff is not None else None

    def calculate(self, atoms=None, properties=('energy', 'forces'), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)

        positions = self.atoms.get_positions()  # shape (N, 3)
        n = len(positions)
        forces = np.zeros((n, 3), dtype=float)
        energy = 0.0

        # Simple O(N^2) loop over unique pairs
        for i in range(n - 1):
            ri = positions[i]
            for j in range(i + 1, n):
                rvec = positions[j] - ri
                r = np.linalg.norm(rvec)
                if r == 0.0:
                    # If two particles coincide, skip to avoid division by zero.
                    # Alternatively, treat as a tiny distance.
                    continue
                if self.cutoff is not None and r > self.cutoff:
                    continue

                dr = r - self.r0
                e_ij = 0.5 * self.k * dr * dr
                energy += e_ij

                f_mag_over_r = self.k * dr / r  # derivative along r̂
                fij = f_mag_over_r * rvec       # force on j due to i
                forces[i] += fij
                forces[j] -= fij                # Newton's 3rd law

        self.results['energy'] = energy
        self.results['forces'] = forces


# System parameters

box_l = 10

# Integration parameters
system = espressomd.System(box_l=[box_l] * 3)

system.time_step = 0.01
system.cell_system.skin = 0.4

p1 = system.part.add(pos=[0,0,0])
p2 = system.part.add(pos=[-0.92,0,0])


## Setup of ase interface and Lennard-Jones calculator

# Mapping of ESPResSo types to ASE types

type_mapping = {0:0} 

ase = espressomd.plugins.ase.ASEInterface(
     system, type_mapping, system.part.all())


harmonic = HarmonicAllPairs(k=1,r0=0)
system.integrator.set_symplectic_euler()
system.thermostat.set_langevin(kT=1,gamma=1,seed=1)
E_kins=[]
E_pots=[]
ase.integrate(100,harmonic)
for i in range(100000):
   ase.integrate(10,harmonic)
   E_pot = 1/2 *system.distance(p1,p2)**2
   E_kin = 1/2 *np.sum(p1.v**2 +p2.v**2)
   E_kins.append(E_kin)
   E_pots.append(E_pot)
   if i%100==0: print(np.mean(E_kins),np.mean(E_pots))

