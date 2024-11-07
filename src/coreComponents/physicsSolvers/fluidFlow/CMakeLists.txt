# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     fluidFlow/CompositionalMultiphaseBase.hpp
     fluidFlow/CompositionalMultiphaseBaseFields.hpp
     fluidFlow/CompositionalMultiphaseStatistics.hpp
     fluidFlow/kernels/IsothermalCompositionalMultiphaseBaseKernels.hpp
     fluidFlow/kernels/ThermalCompositionalMultiphaseBaseKernels.hpp
     fluidFlow/CompositionalMultiphaseFVM.hpp
     fluidFlow/kernels/IsothermalCompositionalMultiphaseFVMKernels.hpp
     fluidFlow/kernels/IsothermalCompositionalMultiphaseFVMKernelUtilities.hpp
     fluidFlow/kernels/ThermalCompositionalMultiphaseFVMKernels.hpp
     fluidFlow/CompositionalMultiphaseHybridFVM.hpp
     fluidFlow/kernels/CompositionalMultiphaseHybridFVMKernels.hpp
     fluidFlow/CompositionalMultiphaseUtilities.hpp
     fluidFlow/ReactiveCompositionalMultiphaseOBL.hpp
     fluidFlow/ReactiveCompositionalMultiphaseOBLFields.hpp
     fluidFlow/kernels/ReactiveCompositionalMultiphaseOBLKernels.hpp
     fluidFlow/FlowSolverBase.hpp
     fluidFlow/FlowSolverBaseFields.hpp
     fluidFlow/kernels/FlowSolverBaseKernels.hpp
     fluidFlow/kernels/FluxKernelsHelper.hpp
     fluidFlow/kernels/HybridFVMHelperKernels.hpp
     fluidFlow/SourceFluxStatistics.hpp
     fluidFlow/proppantTransport/ProppantTransport.hpp
     fluidFlow/proppantTransport/ProppantTransportFields.hpp
     fluidFlow/proppantTransport/ProppantTransportKernels.hpp
     fluidFlow/SinglePhaseBase.hpp
     fluidFlow/SinglePhaseBaseFields.hpp
     fluidFlow/kernels/SinglePhaseBaseKernels.hpp
     fluidFlow/SinglePhaseStatistics.hpp
     fluidFlow/SinglePhaseFVM.hpp
     fluidFlow/kernels/SinglePhaseFVMKernels.hpp
     fluidFlow/SinglePhaseHybridFVM.hpp
     fluidFlow/kernels/SinglePhaseHybridFVMKernels.hpp
     fluidFlow/SinglePhaseProppantBase.hpp
     fluidFlow/kernels/SinglePhaseProppantBaseKernels.hpp
     fluidFlow/kernels/SinglePhaseProppantFluxKernels.hpp
     fluidFlow/kernels/StabilizedCompositionalMultiphaseFVMKernels.hpp
     fluidFlow/kernels/StabilizedSinglePhaseFVMKernels.hpp
     fluidFlow/StencilAccessors.hpp
     fluidFlow/StencilDataCollection.hpp
     fluidFlow/kernels/ThermalSinglePhaseBaseKernels.hpp
     fluidFlow/kernels/ThermalSinglePhaseFVMKernels.hpp
     fluidFlow/LogLevelsInfo.hpp
     fluidFlow/wells/CompositionalMultiphaseWell.hpp
     fluidFlow/wells/CompositionalMultiphaseWellFields.hpp
     fluidFlow/wells/kernels/CompositionalMultiphaseWellKernels.hpp
     fluidFlow/wells/SinglePhaseWell.hpp
     fluidFlow/wells/SinglePhaseWellFields.hpp
     fluidFlow/wells/kernels/SinglePhaseWellKernels.hpp
     fluidFlow/wells/WellConstants.hpp
     fluidFlow/wells/WellControls.hpp
     fluidFlow/wells/WellSolverBase.hpp
     fluidFlow/wells/WellSolverBaseFields.hpp
     fluidFlow/wells/LogLevelsInfo.hpp
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     fluidFlow/CompositionalMultiphaseBase.cpp
     fluidFlow/CompositionalMultiphaseFVM.cpp
     fluidFlow/CompositionalMultiphaseStatistics.cpp
     fluidFlow/kernels/IsothermalCompositionalMultiphaseFVMKernels.cpp
     fluidFlow/CompositionalMultiphaseHybridFVM.cpp
     fluidFlow/kernels/CompositionalMultiphaseHybridFVMKernels.cpp
     fluidFlow/ReactiveCompositionalMultiphaseOBL.cpp
     fluidFlow/FlowSolverBase.cpp
     fluidFlow/proppantTransport/ProppantTransport.cpp
     fluidFlow/proppantTransport/ProppantTransportKernels.cpp
     fluidFlow/SinglePhaseBase.cpp
     fluidFlow/SinglePhaseStatistics.cpp
     fluidFlow/SinglePhaseFVM.cpp
     fluidFlow/SinglePhaseHybridFVM.cpp
     fluidFlow/SinglePhaseProppantBase.cpp
     fluidFlow/kernels/SinglePhaseProppantFluxKernels.cpp
     fluidFlow/SourceFluxStatistics.cpp
     fluidFlow/StencilDataCollection.cpp
     fluidFlow/wells/CompositionalMultiphaseWell.cpp
     fluidFlow/wells/kernels/CompositionalMultiphaseWellKernels.cpp
     fluidFlow/wells/SinglePhaseWell.cpp
     fluidFlow/wells/kernels/SinglePhaseWellKernels.cpp
     fluidFlow/wells/WellControls.cpp
     fluidFlow/wells/WellSolverBase.cpp 
     PARENT_SCOPE )


