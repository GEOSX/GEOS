# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 TotalEnergies
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
#
#--------------------------------------------------------------------------------------------------

#[[
Package: events

Contains GEOS event types.
Manages the events and tasks.
#]]

#
# Specify all headers
#
set( events_headers
     EventBase.hpp
     EventManager.hpp
     HaltEvent.hpp
     PeriodicEvent.hpp
     SoloEvent.hpp
     tasks/TaskBase.hpp
     tasks/TasksManager.hpp
     mpmEvents/MPMEventManager.hpp
     mpmEvents/MPMEventBase.hpp
     mpmEvents/MaterialSwapMPMEvent.hpp
     mpmEvents/AnnealMPMEvent.hpp
     mpmEvents/HealMPMEvent.hpp
     mpmEvents/CrystalHealMPMEvent.hpp
     mpmEvents/InsertPeriodicContactSurfacesMPMEvent.hpp
     mpmEvents/MachineSampleMPMEvent.hpp
     mpmEvents/FrictionCoefficientSwapMPMEvent.hpp
     mpmEvents/BodyForceUpdateMPMEvent.hpp
     mpmEvents/DeformationUpdateMPMEvent.hpp
     mpmEvents/BoreholePressureMPMEvent.hpp
     mpmEvents/CohesiveZoneReferenceMPMEvent.hpp
   )

#
# Specify all sources
#
set( events_sources
     EventBase.cpp
     EventManager.cpp
     HaltEvent.cpp
     PeriodicEvent.cpp
     SoloEvent.cpp
     tasks/TaskBase.cpp
     tasks/TasksManager.cpp
     mpmEvents/MPMEventManager.cpp
     mpmEvents/MPMEventBase.cpp
     mpmEvents/MaterialSwapMPMEvent.cpp
     mpmEvents/AnnealMPMEvent.cpp
     mpmEvents/HealMPMEvent.cpp
     mpmEvents/CrystalHealMPMEvent.cpp
     mpmEvents/InsertPeriodicContactSurfacesMPMEvent.cpp
     mpmEvents/MachineSampleMPMEvent.cpp
     mpmEvents/FrictionCoefficientSwapMPMEvent.cpp
     mpmEvents/BodyForceUpdateMPMEvent.cpp
     mpmEvents/DeformationUpdateMPMEvent.cpp
     mpmEvents/BoreholePressureMPMEvent.cpp
     mpmEvents/CohesiveZoneReferenceMPMEvent.cpp
   )

set( dependencyList ${parallelDeps} functions )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       events
                 SOURCES    ${events_sources}
                 HEADERS    ${events_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( events PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS events LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
