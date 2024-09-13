#
# Specify all headers
#
set( mainInterface_headers
     GeosxState.hpp
     ProblemManager.hpp
     initialization.hpp
     version.hpp
   )

#
# Specify all sources
#
set( mainInterface_sources
     GeosxState.cpp
     ProblemManager.cpp
     initialization.cpp
     version.cpp
   )

set( dependencyList ${parallelDeps} physicsSolvers constitutiveDrivers optionparser )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       mainInterface
                 SOURCES    ${mainInterface_sources}
                 HEADERS    ${mainInterface_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

add_dependencies( mainInterface generate_version )

target_include_directories( mainInterface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS mainInterface LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
