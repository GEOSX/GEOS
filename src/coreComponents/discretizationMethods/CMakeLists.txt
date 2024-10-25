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

geos_decorate_link_dependencies( LIST decoratedDependencies
     DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       discretizationMethods
                 SOURCES    ${mainInterface_sources}
                 HEADERS    ${mainInterface_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( discretizationMethods PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS discretizationMethods LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
