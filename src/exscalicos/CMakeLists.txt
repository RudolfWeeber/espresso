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

# The library is compiled into the ESPResSo core, since Kokkos, Cabana and
# heFFTe are static libraries that must not appear in multiple shared objects
add_library(espresso_exscalicos OBJECT src/Coulomb.cpp)
add_library(espresso::exscalicos ALIAS espresso_exscalicos)
espresso_set_common_target_properties(espresso_exscalicos)
set_target_properties(espresso_exscalicos PROPERTIES POSITION_INDEPENDENT_CODE
                                                     ON)
target_link_libraries(
  espresso_exscalicos PUBLIC MPI::MPI_CXX
  PRIVATE Exscalicos::Exscalicos Kokkos::kokkos Cabana::Core Cabana::Grid
          Heffte::Heffte espresso::compiler_flags)
espresso_set_kokkos_properties_on_source_files(espresso_exscalicos)

# The bridge header uses the short-range kernels of the library, which only
# depend on Kokkos and are evaluated in the ESPResSo near-field loop
target_include_directories(
  espresso_exscalicos
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<TARGET_PROPERTY:Exscalicos::Exscalicos,INTERFACE_INCLUDE_DIRECTORIES>
)
