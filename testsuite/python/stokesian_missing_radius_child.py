#
# Copyright (C) 2026 The ESPResSo project
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
Helper script for ``propagation_stokesian.py``.

Sets up Stokesian Dynamics with a radius for type 0 only, but adds a
particle of type 1 (no radius defined). On the unfixed core this makes
rank 0 throw ``std::out_of_range`` from ``radii.at(p.type)`` *after* the
collective gather but *before* the collective scatter, so the other ranks
block forever in ``MPI_Scatterv`` -> MPI deadlock.

A correct core must raise a coordinated ESPResSo runtime error on all
ranks instead of hanging.

The script is meant to be launched under ``mpiexec`` with a hard timeout.
A clean coordinated error -> exit code 0 ("OK" printed). A hang/timeout
-> the launcher kills it (the deadlock signature).
"""
import espressomd

system = espressomd.System(box_l=[10.0, 10.0, 10.0])
system.cell_system.skin = 0.0
system.periodicity = [False, False, False]
system.time_step = 0.01
system.min_global_cut = 4.0

# place two particles in different regions so they can land on
# different MPI ranks under the default domain decomposition
system.part.add(pos=[1.0, 1.0, 1.0], type=0)
system.part.add(pos=[9.0, 9.0, 9.0], type=1)

system.thermostat.set_stokesian(kT=0.0, seed=42)
# radius defined for type 0 only; type 1 has no radius entry
system.integrator.set_stokesian_dynamics(viscosity=1.0, radii={0: 1.0})

try:
    system.integrator.run(1)
except Exception as err:
    # a clean, coordinated error is the desired behaviour
    if "radius" in str(err).lower():
        print("OK: coordinated runtime error raised")
    else:
        print(f"UNEXPECTED ERROR: {err}")
        raise
else:
    print("NO ERROR RAISED (integration completed unexpectedly)")
