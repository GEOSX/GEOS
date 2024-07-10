#
# Specify all headers
#
set( testingUtilities_headers
     TestingTasks.hpp
   )

#
# Specify all sources
#
set( testingUtilities_sources
     TestingTasks.cpp
   )

#
# Specify all dependencies
#
set( dependencyList ${parallelDeps} gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  list( APPEND dependencyList geosx_core )
else()
  list( APPEND dependencyList ${geosx_core_libs} )
endif()
   
blt_add_library( NAME       testingUtilities
                 SOURCES    ${testingUtilities_sources}
                 HEADERS    ${testingUtilities_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( testingUtilities PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
