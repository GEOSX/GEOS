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
set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList events )
geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       testingUtilities
                 SOURCES    ${testingUtilities_sources}
                 HEADERS    ${testingUtilities_headers}
                 DEPENDS_ON ${decoratedDependencies} ${tplDependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( testingUtilities PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
target_include_directories( testingUtilities PUBLIC ${CMAKE_BINARY_DIR}/include )
