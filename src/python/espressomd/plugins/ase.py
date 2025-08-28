#
# Copyright (C) 2024 The ESPResSo project
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

import dataclasses
import typing
import ase
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np
import espressomd.code_info
if typing.TYPE_CHECKING:
    from espressomd.system import System
    from espressomd.particle_data import ParticleSlice


class ASEInterface:
    """
    ASE interface for ESPResSo with enhanced functionality for calculator integration.
    """

    def __init__(self, system: "System", type_mapping: dict, 
                 particle_slice: "ParticleSlice", 
                 export_charges: bool = False, export_masses: bool = False, 
                 export_momenta: bool = False):
        """
        Initialize ASE interface.

        Parameters
        ----------
        system : espressomd.system.System
            The ESPResSo system object
        type_mapping : dict
            Mapping of ESPResSo particle types to ASE symbols. E.g. ``{0: "H", 1: "O"}``
        particle_slice : espressomd.particle_data.ParticleSlice
            The particle slice to work on
        export_charges : bool, optional
            Whether to make particle charges available to ASE
        export_masses : bool, optional
            Whether to make particle masses available to ASE
        export_momenta : bool, optional
            Whether to make particle momenta available to ASE
        """
        # Check that EXTERNAL_FORCES feature is available
        if "EXTERNAL_FORCES" not in espressomd.code_info.features():
            raise RuntimeError(
                "ASE interface requires EXTERNAL_FORCES feature")

        self.type_mapping = type_mapping
        self._system = system
        self.particle_slice = particle_slice
        self.export_charges = export_charges
        self.export_masses = export_masses
        self.export_momenta = export_momenta
        self.calculator = None
        self.atoms = None

        self.reset()
        self.update_ase()

    def __getstate__(self):
        return {"type_mapping": self.type_mapping}

    def reset(self):
        """
        Re-create the ASE atoms object using ESPResSo system properties.

        Uses the current particle slice to create a new ASE atoms object with
        positions, types, periodicity and box dimensions.
        """
        particles = self.particle_slice
        positions = np.copy(particles.pos)
        types = np.copy(particles.type)

        # Check that all types are in the type mapping
        unknown_types = set(types) - set(self.type_mapping)
        if unknown_types:
            raise RuntimeError(
                f"Particle types '{
                    unknown_types}' haven't been registered in the ASE type map"
            )

        # Check for virtual sites
        if any(p.is_virtual() for p in particles):
            raise RuntimeError("ASE doesn't support virtual sites")

        # Create new atoms object
        self.atoms = ase.Atoms(
            positions=positions,
            symbols=[self.type_mapping[t] for t in types],
            pbc=np.copy(self._system.periodicity),
            cell=np.copy(self._system.box_l),
        )
        self.update_ase()

    def update_ase(self, skip_charge_update: bool = False, skip_mass_update: bool = False):
          """
          Update the arrays in the atoms object based on the desired properties.
  
          Parameters
          ----------
          skip_charge_update : bool, optional
              Whether to skip updating charges
          skip_mass_update : bool, optional
              Whether to skip updating masses
          """
          if self.atoms is None:
            raise RuntimeError(
              "atoms object not initialized, call reset() first")
  
          particles = self.particle_slice
  
          # Always update positions (pos -> positions)
          self.atoms.positions = np.copy(particles.pos)
  
          # Update charges if requested and not skipped
          if self.export_charges and not skip_charge_update:
              charges = np.copy(particles.q)
              self.atoms.set_initial_charges(charges)
  
          # Update masses if requested and not skipped  
          if self.export_masses and not skip_mass_update:
              masses = np.copy(particles.mass)
              self.atoms.set_masses(masses)
  
          # Update momenta if requested
          if self.export_momenta:
              momenta = np.copy(particles.v) * \
                            np.copy(particles.mass)[:, np.newaxis]
              self.atoms.set_momenta(momenta)

    def integrate(self, steps: int, skip_charge_update: bool = False, 
                 skip_mass_update: bool = False) -> int:
                  """
        Integrate the system for the specified number of steps.

        For each step:
        1. Update ASE with current particle data
        2. Get forces from calculator
        3. Set external forces on particles
        4. Run one integration step

        Parameters
        ----------
        steps : int
            Number of integration steps to perform
        skip_charge_update : bool, optional
            Whether to skip charge updates during ASE updates
        skip_mass_update : bool, optional
            Whether to skip mass updates during ASE updates

        Returns
        -------
        int
            Number of steps actually performed
        """
                  if self.calculator is None:
                  raise RuntimeError(
                      "No calculator assigned. Set self.calculator before integrating.")

                  for step in range(steps):
                  # Update ASE with current particle data
                  self.update_ase(skip_charge_update=skip_charge_update, 
                           skip_mass_update=skip_mass_update)

            # Get forces from ASE calculator
            forces = self.calculator.get_forces(self.atoms)

            # Set external forces on particle slice
            self.particle_slice.ext_force = forces

            # Run one integration step
            self._system.integrator.run(1)

        return steps
