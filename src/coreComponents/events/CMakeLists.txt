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
