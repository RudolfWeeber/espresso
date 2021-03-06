#
# Copyright (C) 2013,2014 The ESPResSo project
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
from utils cimport *


def integrate(nSteps, recalc_forces=False, reuse_forces=False):
  checkTypeOrExcept(nSteps,1,int,"Integrate requires a positive integer")
  checkTypeOrExcept(recalc_forces,1,bool,"Integrate requires a positive integer")
  checkTypeOrExcept(nSteps,1,int,"Integrate requires a positive integer")

  if (python_integrate(nSteps, recalc_forces, reuse_forces)):
    print ( mpiRuntimeErrorCollectorGather() )
    raise Exception("Encoutered errors during integrate")


def setIntegratorNVT():
  integrate_set_nvt()

def setIntegratorIsotropicNPT(ext_pressure=0.0, piston=0.0, xdir=0, ydir=0, zdir=0, cubic_box=False):
  checkTypeOrExcept(ext_pressure,1,float,"NPT parameter ext_pressure must be a float")
  checkTypeOrExcept(piston,1,float,"NPT parameter piston must be a float")
  checkTypeOrExcept(xdir,1,int,"NPT parameter xdir must be an int")
  checkTypeOrExcept(ydir,1,int,"NPT parameter ydir must be an int")
  checkTypeOrExcept(zdir,1,int,"NPT parameter zdir must be an int")
  if (integrate_set_npt_isotropic(ext_pressure, piston, xdir, ydir, zdir, cubic_box)):
    print ( mpiRuntimeErrorCollectorGather() )
    raise Exception("Encoutered errors setting up the NPT integrator")
