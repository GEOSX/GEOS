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

set( dependencyList ${parallelDeps} common fileIO dataRepository )

blt_add_library( NAME       events
                 SOURCES    ${events_sources}
                 HEADERS    ${events_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
               )

target_include_directories( events PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

