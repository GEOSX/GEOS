# Specify solver headers
set( physicsSolvers_headers
     LinearSolverParameters.hpp
     NonlinearSolverParameters.hpp
     PhysicsSolverManager.hpp
     SolverBase.hpp
     SolverBaseKernels.hpp
     SolverStatistics.hpp
     FieldStatisticsBase.hpp
     LogLevelsInfo.hpp )

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
endif()

if( GEOS_ENABLE_SURFACEGENERATION )
  add_subdirectory( surfaceGeneration )
endif()

if( GEOS_ENABLE_WAVEPROPAGATION )
  add_subdirectory( wavePropagation )
endif()


if(NOT DEFINED GEOS_KERNEL_SPEC)
  set(GEOS_KERNEL_SPEC "${CMAKE_CURRENT_LIST_DIR}/kernels-all.json" )
endif()
 
file( READ ${GEOS_KERNEL_SPEC} kernelSpecs )
 
set( kernelTemplateFileList
     solidMechanics/kernels/SolidMechanicsKernels.cpp.template
     solidMechanics/kernels/SolidMechanicsFixedStressThermoPoromechanicsKernels.cpp.template
     multiphysics/poromechanicsKernels/PoromechanicsKernels.cpp.template
     multiphysics/poromechanicsKernels/PoromechanicsEFEMKernels.cpp.template
     multiphysics/poromechanicsKernels/ThermoPoromechanicsKernels.cpp.template
     fluidFlow/SinglePhaseHybridFVMKernels.cpp.template
     fluidFlow/CompositionalMultiphaseHybridFVMKernels_upwinding.cpp.template
     fluidFlow/CompositionalMultiphaseHybridFVMKernels_assembly.cpp.template
     fluidFlow/CompositionalMultiphaseHybridFVMKernels_flux.cpp.template
     fluidFlow/ReactiveCompositionalMultiphaseOBLKernels.cpp.template
  )
 
foreach( kernelTemplateFile ${kernelTemplateFileList} )
  get_filename_component( jsonKey ${kernelTemplateFile} NAME_WE )
  generateKernels( TEMPLATE ${kernelTemplateFile}
                   JSON kernelSpecs
                   KEY ${jsonKey}
                   HEADERS headerFiles
                   SOURCES sourceFiles )

  list(APPEND physicsSolvers_headers ${headerFiles})
  list(APPEND physicsSolvers_sources ${sourceFiles})
endforeach()
 
# TODO: The two kernels below have non-matching file names and JSON keys.
#       Either fix them to follow pattern, or come up with another mechanism.

generateKernels( TEMPLATE "fluidFlow/CompositionalMultiphaseFVMKernel_isothermal.cpp.template"
                JSON kernelSpecs
                KEY "IsothermalCompositionalMultiphaseFVM_DirichletFaceBasedAssemblyKernels"
                HEADERS headerFiles
                SOURCES sourceFiles )
list(APPEND physicsSolvers_headers ${headerFiles})
list(APPEND physicsSolvers_sources ${sourceFiles})

generateKernels( TEMPLATE "fluidFlow/CompositionalMultiphaseFVMKernel_thermal.cpp.template"
                JSON kernelSpecs
                KEY "ThermalCompositionalMultiphaseFVM_DirichletFaceBasedAssemblyKernels"
                HEADERS headerFiles
                SOURCES sourceFiles )
list(APPEND physicsSolvers_headers ${headerFiles})
list(APPEND physicsSolvers_sources ${sourceFiles})


set( dependencyList ${parallelDeps} fileIO discretizationMethods events linearAlgebra )

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


