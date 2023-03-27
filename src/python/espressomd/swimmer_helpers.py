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

import numpy as np
from .math import calc_quaternions_from_angles
from .virtual_sites import VirtualSitesRelative

def add_dipole_particle(system, particle, dipole_length, dipole_particle_type, mode = "pusher"):
    """
    Adds a virtual site pointing opposite of the swimmer particle 
    either in front of (mode="puller") or behind (mode="pusher") the swimmer.
    
    Parameters
    ----------
    
    system : :obj:`espressomd.System`
        The system that the particle belongs to.
    particle : :obj:`espressomd.ParticleHandle`
        The swimmer particle that will be paired with the virtual dipole particle.
        Must have the ``swimming`` attribute already set up to get the correct value for ``f_swim``
    dipole_length : :obj:`float`
        The distance between the swimmer and the virtual dipole particle.
    dipole_particle_type : :obj:`int`
        The type of the virtual dipole particle.
    mode : :obj:`str`
        Allowed values: "pusher" (default) and "puller". Determines wether the virtual dipole particle 
        will be placed in front of or behind the swimmer.
    """
    if mode == "pusher":
        dip_sign = -1
    elif mode == "puller":
        dip_sign = 1
    else:
        raise ValueError(f"'mode' must be 'pusher' or 'puller'. You gave '{mode}'")
    
    if dipole_length < 0:
        raise ValueError("'dipole_length' must be >= 0.")
    
    if not isinstance(system.virtual_sites, VirtualSitesRelative):
        raise RuntimeError("system.virtual_sites must be espressomd.virtual_sites.VirtualSitesRelative.")
    if not system.virtual_sites.have_quaternion :
        raise RuntimeError("system.virtual_sites must have quaternion option turned on ('have_quaternion = True').")

    dip_partcl = system.part.add(pos = particle.pos + dip_sign * dipole_length * particle.director,
                           virtual = True,
                           type=dipole_particle_type,
                            swimming={
                                "f_swim": particle.swimming["f_swim"],
                                "is_engine_force_applier": True,
                            },
        )
    dip_partcl.vs_auto_relate_to(particle)
    dip_partcl.vs_quat = calc_quaternions_from_angles(np.pi, 0)
    return dip_partcl
    