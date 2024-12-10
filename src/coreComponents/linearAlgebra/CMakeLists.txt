# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 Total, S.A
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
#
#--------------------------------------------------------------------------------------------------

#[[
Package: linearAlgebra

Contains:
 - management implementations to handle all degrees of freedom associated with field element on
    mesh elements.
 - linear solvers implementations.
 - interfaces to different external linear solvers libraries (hypre, Petsc, trilinos).
#]]

#
# Specify all headers
#
set( linearAlgebra_headers
     DofManager.hpp
     DofManagerHelpers.hpp
     common/LinearOperator.hpp
     common/LinearSolverBase.hpp
     common/PreconditionerBase.hpp
     common/common.hpp
     common/traits.hpp
     interfaces/InterfaceTypes.hpp
     interfaces/MatrixBase.hpp
     interfaces/VectorBase.hpp
     solvers/BicgstabSolver.hpp
     solvers/BlockPreconditioner.hpp
     solvers/CgSolver.hpp
     solvers/GmresSolver.hpp
     solvers/KrylovSolver.hpp
     solvers/KrylovUtils.hpp
     solvers/PreconditionerBlockJacobi.hpp
     solvers/PreconditionerIdentity.hpp
     solvers/PreconditionerJacobi.hpp
     solvers/SeparateComponentPreconditioner.hpp
     utilities/Arnoldi.hpp
     utilities/BlockOperator.hpp
     utilities/BlockOperatorView.hpp
     utilities/BlockOperatorWrapper.hpp
     utilities/BlockVector.hpp
     utilities/BlockVectorView.hpp
     utilities/BlockVectorWrapper.hpp
     utilities/ComponentMask.hpp
     utilities/InverseNormalOperator.hpp
     utilities/LAIHelperFunctions.hpp
     utilities/LinearSolverParameters.hpp
     utilities/LinearSolverResult.hpp
     utilities/NormalOperator.hpp
     utilities/ReverseCutHillMcKeeOrdering.hpp
     utilities/TransposeOperator.hpp )
#
# Specify all sources
#
set( linearAlgebra_sources
     DofManager.cpp
     solvers/BicgstabSolver.cpp
     solvers/BlockPreconditioner.cpp
     solvers/CgSolver.cpp
     solvers/GmresSolver.cpp
     solvers/KrylovSolver.cpp
     solvers/SeparateComponentPreconditioner.cpp
     utilities/ReverseCutHillMcKeeOrdering.cpp )

set( dependencyList ${parallelDeps} mesh denseLinearAlgebra finiteVolume )
set( tplDependencyList "" )

list( APPEND linearAlgebra_headers interfaces/direct/SuiteSparse.hpp )
list( APPEND linearAlgebra_sources interfaces/direct/SuiteSparse.cpp )
list( APPEND tplDependencyList suitesparse )

if( ENABLE_SUPERLU_DIST )
      list( APPEND linearAlgebra_headers interfaces/direct/SuperLUDist.hpp )
      list( APPEND linearAlgebra_sources interfaces/direct/SuperLUDist.cpp )
      list( APPEND tplDependencyList superlu_dist )
endif( )

if( ENABLE_TRILINOS )

    list( APPEND linearAlgebra_headers
          interfaces/trilinos/EpetraExport.hpp
          interfaces/trilinos/EpetraMatrix.hpp
          interfaces/trilinos/EpetraUtils.hpp
          interfaces/trilinos/EpetraVector.hpp
          interfaces/trilinos/TrilinosInterface.hpp
          interfaces/trilinos/TrilinosPreconditioner.hpp
          interfaces/trilinos/TrilinosSolver.hpp )

    list( APPEND linearAlgebra_sources
          interfaces/trilinos/EpetraExport.cpp
          interfaces/trilinos/EpetraMatrix.cpp
          interfaces/trilinos/EpetraVector.cpp
          interfaces/trilinos/TrilinosInterface.cpp
          interfaces/trilinos/TrilinosPreconditioner.cpp
          interfaces/trilinos/TrilinosSolver.cpp )

    list( APPEND tplDependencyList trilinos )

endif()

if( ENABLE_HYPRE )

    list( APPEND linearAlgebra_headers
          interfaces/hypre/HypreExport.hpp
          interfaces/hypre/HypreInterface.hpp
          interfaces/hypre/HypreKernels.hpp
          interfaces/hypre/HypreMGR.hpp
          interfaces/hypre/HypreMatrix.hpp
          interfaces/hypre/HyprePreconditioner.hpp
          interfaces/hypre/HypreSolver.hpp
          interfaces/hypre/HypreUtils.hpp
          interfaces/hypre/HypreVector.hpp
          interfaces/hypre/mgrStrategies/CompositionalMultiphaseFVM.hpp
          interfaces/hypre/mgrStrategies/CompositionalMultiphaseHybridFVM.hpp
          interfaces/hypre/mgrStrategies/CompositionalMultiphaseReservoirFVM.hpp
          interfaces/hypre/mgrStrategies/CompositionalMultiphaseReservoirHybridFVM.hpp
          interfaces/hypre/mgrStrategies/HybridSinglePhasePoromechanics.hpp
          interfaces/hypre/mgrStrategies/Hydrofracture.hpp
          interfaces/hypre/mgrStrategies/LagrangianContactMechanics.hpp
          interfaces/hypre/mgrStrategies/MultiphasePoromechanics.hpp
          interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsEmbeddedFractures.hpp
          interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsConformingFractures.hpp
          interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsReservoirFVM.hpp
          interfaces/hypre/mgrStrategies/SinglePhaseReservoirFVM.hpp
          interfaces/hypre/mgrStrategies/SinglePhaseReservoirHybridFVM.hpp
          interfaces/hypre/mgrStrategies/SolidMechanicsEmbeddedFractures.hpp
          interfaces/hypre/mgrStrategies/ThermalCompositionalMultiphaseFVM.hpp
          interfaces/hypre/mgrStrategies/ThermalCompositionalMultiphaseReservoirFVM.hpp
          interfaces/hypre/mgrStrategies/ThermalMultiphasePoromechanics.hpp
          interfaces/hypre/mgrStrategies/ThermalSinglePhasePoromechanics.hpp )
    list( APPEND linearAlgebra_sources
          interfaces/hypre/HypreExport.cpp
          interfaces/hypre/HypreInterface.cpp
          interfaces/hypre/HypreKernels.cpp
          interfaces/hypre/HypreMGR.cpp
          interfaces/hypre/HypreMatrix.cpp
          interfaces/hypre/HyprePreconditioner.cpp
          interfaces/hypre/HypreSolver.cpp
          interfaces/hypre/HypreUtils.cpp
          interfaces/hypre/HypreVector.cpp )

    list( APPEND tplDependencyList hypre umpire )

endif()

if( ENABLE_PETSC )

    list( APPEND linearAlgebra_headers
          interfaces/petsc/PetscExport.hpp
          interfaces/petsc/PetscInterface.hpp
          interfaces/petsc/PetscMatrix.hpp
          interfaces/petsc/PetscPreconditioner.hpp
          interfaces/petsc/PetscSolver.hpp
          interfaces/petsc/PetscUtils.hpp
          interfaces/petsc/PetscVector.hpp )

    list( APPEND linearAlgebra_sources
          interfaces/petsc/PetscExport.cpp
          interfaces/petsc/PetscInterface.cpp
          interfaces/petsc/PetscMatrix.cpp
          interfaces/petsc/PetscPreconditioner.cpp
          interfaces/petsc/PetscSolver.cpp
          interfaces/petsc/PetscVector.cpp )

    list( APPEND tplDependencyList petsc )

endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       linearAlgebra
                 SOURCES    ${linearAlgebra_sources}
                 HEADERS    ${linearAlgebra_headers}
                 DEPENDS_ON ${decoratedDependencies} ${tplDependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( linearAlgebra PUBLIC ${CMAKE_CURRENT_LIST_DIR} )

install( TARGETS linearAlgebra LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif( )

