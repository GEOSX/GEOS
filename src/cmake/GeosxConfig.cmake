set( PREPROCESSOR_DEFINES ARRAY_BOUNDS_CHECK
                          CALIPER
                          CHAI
                          CUDA
                          CUDA_NVTOOLSEXT
                          HIP
			  FMT_CONST_FORMATTER_WORKAROUND
                          FORTRAN_MANGLE_NO_UNDERSCORE
                          FPE
                          HYPRE
                          MATHPRESSO
                          METIS
                          MKL
                          MPI
                          PARMETIS
                          PETSC
                          PVTPackage
                          PYGEOSX
                          RAJA
                          SCOTCH
                          SEPARATION_COEFFICIENT
                          SUITESPARSE
                          SUPERLU_DIST
                          TIMERS
                          TOTALVIEW_OUTPUT
                          TRILINOS
                          VTK
                          ${externalComponentsList} )

foreach( DEP in ${PREPROCESSOR_DEFINES} )
    if( ${DEP}_FOUND OR ENABLE_${DEP} OR GEOS_ENABLE_${DEP} )
        set( USE_${DEP} TRUE )
        set( GEOS_USE_${DEP} TRUE )
        set( GEOS_USE_${DEP} TRUE )
	message(STATUS "GEOS_USE_${DEP} = ${GEOS_USE_${DEP}}")
    endif()
endforeach()

set( STRICT_PPD OPENMP )

# only activate these options if they are ENABLED AND FOUND, not if either
foreach( DEP in ${STRICT_PPD} )
    if( ${DEP}_FOUND AND ( ENABLE_${DEP} OR GEOS_ENABLE_${DEP} ) )
        set( USE_${DEP} TRUE )
        set( GEOS_USE_${DEP} TRUE )
        set( GEOS_USE_${DEP} TRUE )
	message(STATUS "GEOS_USE_${DEP} = ${GEOS_USE_${DEP}}")
    endif()
endforeach( )

set( GEOS_USE_HYPRE_DEVICE "GEOS_USE_HYPRE_${ENABLE_HYPRE_DEVICE}" )
message( STATUS "GEOS_USE_HYPRE_DEVICE = ${GEOS_USE_HYPRE_DEVICE}")

set( GEOS_CMAKE_BUILD_TYPE "\"${CMAKE_BUILD_TYPE}\"" )

configure_file( ${CMAKE_SOURCE_DIR}/coreComponents/common/GeosxConfig.hpp.in
                ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp )

install( FILES ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp
         DESTINATION ${CMAKE_INSTALL_PREFIX}/include/common )
