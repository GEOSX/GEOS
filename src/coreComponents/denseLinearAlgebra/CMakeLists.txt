# Specify all headers
set( denseLinearAlgebra_headers
     common/layouts.hpp
     interfaces/blaslapack/BlasLapackFunctions.h
     interfaces/blaslapack/BlasLapackLA.hpp )

# Specify all sources
set( denseLinearAlgebra_sources
     interfaces/blaslapack/BlasLapackLA.cpp )

set( dependencyList codingUtilities blas lapack ${parallelDeps} )

blt_add_library( NAME       denseLinearAlgebra
                 SOURCES    ${denseLinearAlgebra_sources}
                 HEADERS    ${denseLinearAlgebra_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS} )

target_include_directories( denseLinearAlgebra PUBLIC ${CMAKE_CURRENT_LIST_DIR} )

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif( )

