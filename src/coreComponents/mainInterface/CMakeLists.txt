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

set( dependencyList ${parallelDeps} physicsSolvers discretizationMethods fieldSpecification linearAlgebra dataRepository events fileIO optionparser )

blt_add_library( NAME       mainInterface
                 SOURCES    ${mainInterface_sources}
                 HEADERS    ${mainInterface_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
               )

add_dependencies( mainInterface generate_version )

target_include_directories( mainInterface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

