# Specify all headers
set( denseLinearAlgebra_headers
     common/layouts.hpp
     denseLASolvers.hpp
     interfaces/blaslapack/BlasLapackFunctions.h
     interfaces/blaslapack/BlasLapackLA.hpp )

# Specify all sources
set( denseLinearAlgebra_sources
     interfaces/blaslapack/BlasLapackLA.cpp )

set( dependencyList codingUtilities blas lapack ${parallelDeps} )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       denseLinearAlgebra
                 SOURCES    ${denseLinearAlgebra_sources}
                 HEADERS    ${denseLinearAlgebra_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS} 
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( denseLinearAlgebra PUBLIC ${CMAKE_CURRENT_LIST_DIR} )

install( TARGETS denseLinearAlgebra LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif( )

