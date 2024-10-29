#
# Specify all headers
#
set( mainInterface_headers
     NumericalMethodsManager.hpp )

set( dependencyList ${parallelDeps} finiteVolume )

#
# Specify all sources
#

set( mainInterface_sources
     NumericalMethodsManager.cpp )

blt_add_library( NAME       discretizationMethods
                 SOURCES    ${mainInterface_sources}
                 HEADERS    ${mainInterface_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
               )

target_include_directories( discretizationMethods PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

