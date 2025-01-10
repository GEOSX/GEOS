# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     multiphysics/CompositionalMultiphaseReservoirAndWells.hpp
     multiphysics/CoupledReservoirAndWellsBase.hpp
     multiphysics/CoupledSolver.hpp
     multiphysics/PoromechanicsSolver.hpp
     multiphysics/FlowProppantTransportSolver.hpp
     multiphysics/HydrofractureSolver.hpp
     multiphysics/HydrofractureSolverKernels.hpp
     multiphysics/MultiphasePoromechanics.hpp
     multiphysics/PhaseFieldFractureSolver.hpp
     multiphysics/PoromechanicsInitialization.hpp
     multiphysics/PoromechanicsFields.hpp
     multiphysics/PoromechanicsInitialization.hpp
     multiphysics/LogLevelsInfo.hpp
     multiphysics/poromechanicsKernels/MultiphasePoromechanics.hpp
     multiphysics/poromechanicsKernels/MultiphasePoromechanics_impl.hpp
     multiphysics/poromechanicsKernels/PoromechanicsBase.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanics_impl.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanicsConformingFractures.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanicsEFEM.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanicsEFEM_impl.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanicsFractures.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanicsEmbeddedFractures.hpp
     multiphysics/poromechanicsKernels/ThermalMultiphasePoromechanics.hpp
     multiphysics/poromechanicsKernels/ThermalMultiphasePoromechanics_impl.hpp
     multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanics.hpp
     multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanics_impl.hpp
     multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanicsEFEM.hpp
     multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanicsEFEM_impl.hpp
     multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanicsConformingFractures.hpp
     multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanicsEmbeddedFractures.hpp
     multiphysics/SinglePhasePoromechanics.hpp
     multiphysics/SinglePhasePoromechanicsEmbeddedFractures.hpp
     multiphysics/SinglePhasePoromechanicsConformingFractures.hpp
     multiphysics/SinglePhaseReservoirAndWells.hpp
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     multiphysics/CompositionalMultiphaseReservoirAndWells.cpp
     multiphysics/CoupledReservoirAndWellsBase.cpp
     multiphysics/FlowProppantTransportSolver.cpp
     multiphysics/HydrofractureSolver.cpp
     multiphysics/MultiphasePoromechanics.cpp
     multiphysics/PhaseFieldFractureSolver.cpp
     multiphysics/PoromechanicsInitialization.cpp
     multiphysics/SinglePhasePoromechanics.cpp
     multiphysics/SinglePhasePoromechanicsEmbeddedFractures.cpp
     multiphysics/SinglePhasePoromechanicsConformingFractures.cpp
     multiphysics/SinglePhaseReservoirAndWells.cpp
     PARENT_SCOPE )

#include( multiphysics/poromechanicsKernels/PoromechanicsKernels.cmake)
