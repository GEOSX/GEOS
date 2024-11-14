# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 Total, S.A
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
#
#--------------------------------------------------------------------------------------------------

#[[
Package: constitutiveDrivers

Contains fluid and rock constitutive model tests to calibrate constitutive model parameters to
 experimental data.
#]]

#
# Specify all headers
#
set( constitutiveDrivers_headers
     fluid/multiFluid/PVTDriver.hpp
     fluid/multiFluid/PVTDriverRunTest.hpp
     fluid/multiFluid/reactive/ReactiveFluidDriver.hpp
     relativePermeability/RelpermDriver.hpp
     relativePermeability/RelpermDriverRunTest.hpp
     solid/TriaxialDriver.hpp
   )
#
# Specify all sources
#
set( constitutiveDrivers_sources
     fluid/multiFluid/PVTDriver.cpp
     fluid/multiFluid/blackOil/PVTDriverRunTestDeadOilFluid.cpp
     fluid/multiFluid/CO2Brine/PVTDriverRunTestCO2BrinePhillipsFluid.cpp
     fluid/multiFluid/CO2Brine/PVTDriverRunTestCO2BrinePhillipsThermalFluid.cpp
     fluid/multiFluid/CO2Brine/PVTDriverRunTestCO2BrineEzrokhiFluid.cpp
     fluid/multiFluid/CO2Brine/PVTDriverRunTestCO2BrineEzrokhiThermalFluid.cpp
     fluid/multiFluid/compositional/PVTDriverRunTestCompositionalTwoPhaseConstantViscosity.cpp
     fluid/multiFluid/compositional/PVTDriverRunTestCompositionalTwoPhaseLohrenzBrayClarkViscosity.cpp
     fluid/multiFluid/reactive/ReactiveFluidDriver.cpp
     relativePermeability/RelpermDriver.cpp
     relativePermeability/RelpermDriverBrooksCoreyBakerRunTest.cpp
     relativePermeability/RelpermDriverBrooksCoreyStone2RunTest.cpp
     relativePermeability/RelpermDriverBrooksCoreyRunTest.cpp
     relativePermeability/RelpermDriverVanGenuchtenBakerRunTest.cpp
     relativePermeability/RelpermDriverVanGenuchtenStone2RunTest.cpp
     relativePermeability/RelpermDriverTableRelativeRunTest.cpp
     relativePermeability/RelpermDriverTableRelativeHysteresisRunTest.cpp
     solid/TriaxialDriver.cpp
   )

set( dependencyList ${parallelDeps} constitutive events )

if( ENABLE_PVTPackage )
    set( constitutiveDrivers_sources
         ${constitutiveDrivers_sources}
          fluid/multiFluid/compositional/PVTDriverRunTestCompositionalMultiphaseFluid.cpp
          )
endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       constitutiveDrivers
                 SOURCES    ${constitutiveDrivers_sources}
                 HEADERS    ${constitutiveDrivers_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( constitutiveDrivers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
install( TARGETS constitutiveDrivers LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
