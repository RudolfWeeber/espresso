#
# Copyright (C) 2021-2025 The ESPResSo project
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

set(backends ".cpp")
if(WALBERLA_BUILD_WITH_CUDA)
  list(APPEND backends "CUDA.cu")
endif()

function(configure_codegen_kernel kernel)
  foreach(precision Single Double)
    foreach(backend ${backends})
      set(source_file ${kernel}${precision}Precision${backend})
      target_sources(espresso_walberla_codegen PRIVATE ${source_file})
      espresso_set_walberla_codegen_kernel_source_properties(${source_file})
    endforeach()
  endforeach()
endfunction()

configure_codegen_kernel(ColorGradientInitialPDFsSetter)
if(ESPRESSO_BUILD_WITH_WALBERLA_AVX)
  list(REMOVE_ITEM backends ".cpp")
  list(APPEND backends "AVX.cpp")
endif()
configure_codegen_kernel(ColorGradientStreamSweep)
configure_codegen_kernel(ColorGradientCollideSweep)