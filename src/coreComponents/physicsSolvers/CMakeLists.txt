# Specify solver headers
set( physicsSolvers_headers
     LinearSolverParameters.hpp
     NonlinearSolverParameters.hpp
     PhysicsSolverManager.hpp
     SolverBase.hpp
     SolverBaseKernels.hpp
     SolverStatistics.hpp
     FieldStatisticsBase.hpp
 )

# Specify solver sources
set( physicsSolvers_sources
     LinearSolverParameters.cpp
     NonlinearSolverParameters.cpp
     PhysicsSolverManager.cpp
     SolverBase.cpp
     SolverStatistics.cpp )

if( GEOS_ENABLE_CONTACT )
  add_subdirectory( contact )
endif()

if( GEOS_ENABLE_FLUIDFLOW )
  add_subdirectory( fluidFlow )
endif()

if( GEOS_ENABLE_INDUCEDSEISMICITY )
  add_subdirectory( inducedSeismicity )
endif()

if( GEOS_ENABLE_MULTIPHYSICS )
  add_subdirectory( multiphysics )
  include( multiphysics/poromechanicsKernels/PoromechanicsKernels.cmake)
endif()

if( GEOS_ENABLE_SIMPLEPDE )
  add_subdirectory( simplePDE )
endif()
  
if( GEOS_ENABLE_SOLIDMECHANICS )
  add_subdirectory( solidMechanics )
  include( solidMechanics/kernels/SolidMechanicsKernels.cmake)
endif()

if( GEOS_ENABLE_SURFACEGENERATION )
  add_subdirectory( surfaceGeneration )
endif()

if( GEOS_ENABLE_WAVEPROPAGATION )
  add_subdirectory( wavePropagation )
endif()

set( dependencyList ${parallelDeps} fileIO discretizationMethods events )

if( ENABLE_PYGEOSX )
  list( APPEND physicsSolvers_headers
  python/PySolver.hpp
  python/PySolverType.hpp )
  list( APPEND physicsSolvers_sources
        python/PySolver.cpp )
  list( APPEND dependencyList Python3::Python pylvarray )
endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       physicsSolvers
                 SOURCES    ${physicsSolvers_sources}
                 HEADERS    ${physicsSolvers_headers}
                 DEPENDS_ON ${decoratedDependencies} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS physicsSolvers LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( externalComponentDeps )
  target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()


