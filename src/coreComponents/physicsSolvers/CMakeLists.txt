# Specify solver headers
set( physicsSolvers_headers
     LinearSolverParameters.hpp
     NonlinearSolverParameters.hpp
     PhysicsSolverManager.hpp
     SolverBase.hpp
     SolverBaseKernels.hpp
     SolverStatistics.hpp
     FieldStatisticsBase.hpp
     contact/ContactSolverBase.hpp
     contact/ContactFields.hpp
     contact/SolidMechanicsEFEMKernelsBase.hpp
     contact/SolidMechanicsEFEMKernels.hpp
     contact/SolidMechanicsEFEMStaticCondensationKernels.hpp
     contact/SolidMechanicsEFEMKernelsHelper.hpp
     contact/SolidMechanicsEmbeddedFractures.hpp
     contact/SolidMechanicsLagrangeContact.hpp
     contact/SolidMechanicsAugmentedLagrangianContact.hpp
     contact/SolidMechanicsALMKernelsBase.hpp
     contact/SolidMechanicsALMKernels.hpp
     contact/SolidMechanicsALMKernelsHelper.hpp
     contact/SolidMechanicsALMJumpUpdateKernels.hpp
     contact/SolidMechanicsALMBubbleKernels.hpp
     fluidFlow/CompositionalMultiphaseBase.hpp
     fluidFlow/CompositionalMultiphaseBaseFields.hpp
     fluidFlow/CompositionalMultiphaseStatistics.hpp
     fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp
     fluidFlow/ThermalCompositionalMultiphaseBaseKernels.hpp
     fluidFlow/CompositionalMultiphaseFVM.hpp
     fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp
     fluidFlow/IsothermalCompositionalMultiphaseFVMKernelUtilities.hpp
     fluidFlow/ThermalCompositionalMultiphaseFVMKernels.hpp
     fluidFlow/CompositionalMultiphaseHybridFVM.hpp
     fluidFlow/CompositionalMultiphaseHybridFVMKernels.hpp
     fluidFlow/CompositionalMultiphaseUtilities.hpp
     fluidFlow/ReactiveCompositionalMultiphaseOBL.hpp
     fluidFlow/ReactiveCompositionalMultiphaseOBLFields.hpp
     fluidFlow/ReactiveCompositionalMultiphaseOBLKernels.hpp
     fluidFlow/FlowSolverBase.hpp
     fluidFlow/FlowSolverBaseFields.hpp
     fluidFlow/FlowSolverBaseKernels.hpp
     fluidFlow/FluxKernelsHelper.hpp
     fluidFlow/HybridFVMHelperKernels.hpp
     fluidFlow/proppantTransport/ProppantTransport.hpp
     fluidFlow/proppantTransport/ProppantTransportFields.hpp
     fluidFlow/proppantTransport/ProppantTransportKernels.hpp
     fluidFlow/SinglePhaseBase.hpp
     fluidFlow/SinglePhaseBaseFields.hpp
     fluidFlow/SinglePhaseBaseKernels.hpp
     fluidFlow/SinglePhaseStatistics.hpp
     fluidFlow/SinglePhaseFVM.hpp
     fluidFlow/SinglePhaseFVMKernels.hpp
     fluidFlow/SinglePhaseHybridFVM.hpp
     fluidFlow/SinglePhaseHybridFVMKernels.hpp
     fluidFlow/SinglePhaseProppantBase.hpp
     fluidFlow/SinglePhaseProppantBaseKernels.hpp
     fluidFlow/SinglePhaseProppantFluxKernels.hpp
     fluidFlow/StabilizedCompositionalMultiphaseFVMKernels.hpp
     fluidFlow/StabilizedSinglePhaseFVMKernels.hpp
     fluidFlow/StencilAccessors.hpp
     fluidFlow/ThermalSinglePhaseBaseKernels.hpp
     fluidFlow/ThermalSinglePhaseFVMKernels.hpp
     fluidFlow/wells/CompositionalMultiphaseWell.hpp
     fluidFlow/wells/CompositionalMultiphaseWellFields.hpp
     fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp
     fluidFlow/wells/SinglePhaseWell.hpp
     fluidFlow/wells/SinglePhaseWellFields.hpp
     fluidFlow/wells/SinglePhaseWellKernels.hpp
     fluidFlow/wells/WellConstants.hpp
     fluidFlow/wells/WellControls.hpp
     fluidFlow/wells/WellSolverBase.hpp
     fluidFlow/wells/WellSolverBaseFields.hpp
     inducedSeismicity/inducedSeismicityFields.hpp
     inducedSeismicity/SeismicityRate.hpp
     inducedSeismicity/SeismicityRateKernels.hpp
     multiphysics/CompositionalMultiphaseReservoirAndWells.hpp
     multiphysics/CoupledReservoirAndWellsBase.hpp
     multiphysics/CoupledSolver.hpp
     multiphysics/PoromechanicsSolver.hpp
     multiphysics/FlowProppantTransportSolver.hpp
     multiphysics/HydrofractureSolver.hpp
     multiphysics/HydrofractureSolverKernels.hpp
     multiphysics/MultiphasePoromechanics.hpp
     multiphysics/PhaseFieldFractureSolver.hpp
     multiphysics/PhaseFieldPoromechanicsSolver.hpp
     multiphysics/PoromechanicsInitialization.hpp
     multiphysics/PoromechanicsFields.hpp
     multiphysics/PoromechanicsInitialization.hpp
     multiphysics/poromechanicsKernels/MultiphasePoromechanics.hpp
     multiphysics/poromechanicsKernels/MultiphasePoromechanics_impl.hpp
     multiphysics/poromechanicsKernels/PoromechanicsBase.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanics_impl.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanicsDamage.hpp  
     multiphysics/poromechanicsKernels/SinglePhasePoromechanicsDamage_impl.hpp     
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
     simplePDE/LaplaceBaseH1.hpp
     simplePDE/LaplaceFEM.hpp
     simplePDE/LaplaceFEMKernels.hpp
     simplePDE/PhaseFieldDamageFEM.hpp
     simplePDE/PhaseFieldDamageFEMKernels.hpp
     simplePDE/PhaseFieldPressurizedDamageFEMKernels.hpp
     solidMechanics/SolidMechanicsFields.hpp
     solidMechanics/SolidMechanicsLagrangianFEM.hpp
     solidMechanics/SolidMechanicsLagrangianFEM.hpp
     solidMechanics/SolidMechanicsLagrangianSSLE.hpp
     solidMechanics/kernels/SolidMechanicsLagrangianFEMKernels.hpp
     solidMechanics/SolidMechanicsMPM.hpp
     solidMechanics/MPMSolverFields.hpp
     solidMechanics/kernels/ExplicitFiniteStrain.hpp
     solidMechanics/kernels/ExplicitFiniteStrain_impl.hpp
     solidMechanics/kernels/ExplicitMPM.hpp
     solidMechanics/kernels/ExplicitSmallStrain.hpp
     solidMechanics/kernels/ExplicitSmallStrain_impl.hpp
     solidMechanics/kernels/FixedStressThermoPoromechanics.hpp
     solidMechanics/kernels/FixedStressThermoPoromechanics_impl.hpp
     solidMechanics/kernels/ImplicitSmallStrainNewmark.hpp
     solidMechanics/kernels/ImplicitSmallStrainNewmark_impl.hpp
     solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp
     solidMechanics/kernels/ImplicitSmallStrainQuasiStatic_impl.hpp
     solidMechanics/SolidMechanicsStateReset.hpp
     solidMechanics/SolidMechanicsStatistics.hpp
     surfaceGeneration/EmbeddedSurfaceGenerator.hpp
     surfaceGeneration/EmbeddedSurfacesParallelSynchronization.hpp
     surfaceGeneration/ParallelTopologyChange.hpp
     surfaceGeneration/SurfaceGenerator.hpp
     surfaceGeneration/SurfaceGeneratorFields.hpp
     surfaceGeneration/kernels/surfaceGenerationKernels.hpp
     surfaceGeneration/kernels/surfaceGenerationKernelsHelpers.hpp
     wavePropagation/shared/WaveSolverBase.hpp
     wavePropagation/shared/WaveSolverUtils.hpp
     wavePropagation/shared/PrecomputeSourcesAndReceiversKernel.hpp
     wavePropagation/sem/acoustic/shared/AcousticFields.hpp
     wavePropagation/sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEM.hpp
     wavePropagation/sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEMKernel.hpp
     wavePropagation/sem/acoustic/shared/AcousticTimeSchemeSEMKernel.hpp
     wavePropagation/sem/acoustic/shared/AcousticMatricesSEMKernel.hpp
     wavePropagation/sem/acoustic/secondOrderEqn/isotropic/AcousticPMLSEMKernel.hpp
     wavePropagation/sem/elastic/shared/ElasticFields.hpp
     wavePropagation/sem/elastic/secondOrderEqn/isotropic/ElasticWaveEquationSEM.hpp
     wavePropagation/sem/elastic/secondOrderEqn/isotropic/ElasticWaveEquationSEMKernel.hpp
     wavePropagation/sem/elastic/shared/ElasticTimeSchemeSEMKernel.hpp
     wavePropagation/sem/elastic/firstOrderEqn/isotropic/ElasticFirstOrderWaveEquationSEM.hpp
     wavePropagation/sem/elastic/firstOrderEqn/isotropic/ElasticFirstOrderWaveEquationSEMKernel.hpp
     wavePropagation/sem/acoustic/firstOrderEqn/isotropic/AcousticFirstOrderWaveEquationSEM.hpp
     wavePropagation/sem/acoustic/firstOrderEqn/isotropic/AcousticFirstOrderWaveEquationSEMKernel.hpp
     wavePropagation/sem/elastic/shared/ElasticMatricesSEMKernel.hpp
     wavePropagation/sem/elastic/secondOrderEqn/anisotropic/ElasticVTIFields.hpp
     wavePropagation/sem/elastic/secondOrderEqn/anisotropic/ElasticVTIWaveEquationSEMKernel.hpp
     wavePropagation/sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIFields.hpp
     wavePropagation/sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIWaveEquationSEM.hpp
     wavePropagation/sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIWaveEquationSEMKernel.hpp
     wavePropagation/sem/acoustoelastic/secondOrderEqn/isotropic/AcoustoElasticFields.hpp
     wavePropagation/sem/acoustoelastic/secondOrderEqn/isotropic/AcoustoElasticTimeSchemeSEMKernel.hpp
     wavePropagation/sem/acoustoelastic/secondOrderEqn/isotropic/AcousticElasticWaveEquationSEM.hpp
     wavePropagation/sem/acoustoelastic/secondOrderEqn/isotropic/AcousticElasticWaveEquationSEMKernel.hpp )

# Specify solver sources
set( physicsSolvers_sources
     LinearSolverParameters.cpp
     NonlinearSolverParameters.cpp
     PhysicsSolverManager.cpp
     SolverBase.cpp
     SolverStatistics.cpp
     contact/ContactSolverBase.cpp
     contact/SolidMechanicsEmbeddedFractures.cpp
     contact/SolidMechanicsLagrangeContact.cpp
     contact/SolidMechanicsAugmentedLagrangianContact.cpp
     fluidFlow/CompositionalMultiphaseBase.cpp
     fluidFlow/CompositionalMultiphaseFVM.cpp
     fluidFlow/CompositionalMultiphaseStatistics.cpp
     fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.cpp
     fluidFlow/CompositionalMultiphaseHybridFVM.cpp
     fluidFlow/CompositionalMultiphaseHybridFVMKernels.cpp
     fluidFlow/ReactiveCompositionalMultiphaseOBL.cpp
     fluidFlow/FlowSolverBase.cpp
     fluidFlow/proppantTransport/ProppantTransport.cpp
     fluidFlow/proppantTransport/ProppantTransportKernels.cpp
     fluidFlow/SinglePhaseBase.cpp
     fluidFlow/SinglePhaseStatistics.cpp
     fluidFlow/SinglePhaseFVM.cpp
     fluidFlow/SinglePhaseHybridFVM.cpp
     fluidFlow/SinglePhaseProppantBase.cpp
     fluidFlow/SinglePhaseProppantFluxKernels.cpp
     fluidFlow/wells/CompositionalMultiphaseWell.cpp
     fluidFlow/wells/CompositionalMultiphaseWellKernels.cpp
     fluidFlow/wells/SinglePhaseWell.cpp
     fluidFlow/wells/SinglePhaseWellKernels.cpp
     fluidFlow/wells/WellControls.cpp
     fluidFlow/wells/WellSolverBase.cpp
     inducedSeismicity/SeismicityRate.cpp
     multiphysics/CompositionalMultiphaseReservoirAndWells.cpp
     multiphysics/CoupledReservoirAndWellsBase.cpp
     multiphysics/FlowProppantTransportSolver.cpp
     multiphysics/HydrofractureSolver.cpp
     multiphysics/MultiphasePoromechanics.cpp
     multiphysics/PhaseFieldFractureSolver.cpp
     multiphysics/PhaseFieldPoromechanicsSolver.cpp
     multiphysics/PoromechanicsInitialization.cpp
     multiphysics/SinglePhasePoromechanics.cpp
     multiphysics/SinglePhasePoromechanicsEmbeddedFractures.cpp
     multiphysics/SinglePhasePoromechanicsConformingFractures.cpp
     multiphysics/SinglePhaseReservoirAndWells.cpp
     simplePDE/LaplaceBaseH1.cpp
     simplePDE/LaplaceFEM.cpp
     simplePDE/PhaseFieldDamageFEM.cpp
     solidMechanics/SolidMechanicsLagrangianFEM.cpp
     solidMechanics/SolidMechanicsLagrangianSSLE.cpp
     solidMechanics/SolidMechanicsMPM.cpp
     solidMechanics/SolidMechanicsStateReset.cpp
     solidMechanics/SolidMechanicsStatistics.cpp
     surfaceGeneration/EmbeddedSurfaceGenerator.cpp
     surfaceGeneration/EmbeddedSurfacesParallelSynchronization.cpp
     surfaceGeneration/ParallelTopologyChange.cpp
     surfaceGeneration/SurfaceGenerator.cpp
     wavePropagation/shared/WaveSolverBase.cpp
     wavePropagation/sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEM.cpp
     wavePropagation/sem/elastic/secondOrderEqn/isotropic/ElasticWaveEquationSEM.cpp
     wavePropagation/sem/elastic/firstOrderEqn/isotropic/ElasticFirstOrderWaveEquationSEM.cpp
     wavePropagation/sem/acoustic/firstOrderEqn/isotropic/AcousticFirstOrderWaveEquationSEM.cpp
     wavePropagation/sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIWaveEquationSEM.cpp
     wavePropagation/sem/acoustoelastic/secondOrderEqn/isotropic/AcousticElasticWaveEquationSEM.cpp )

  include( solidMechanics/kernels/SolidMechanicsKernels.cmake)
  include( multiphysics/poromechanicsKernels/PoromechanicsKernels.cmake)

set( dependencyList ${parallelDeps} constitutive mesh linearAlgebra discretizationMethods events )
if( ENABLE_PYGEOSX )
  list( APPEND physicsSolvers_headers
  python/PySolver.hpp
  python/PySolverType.hpp )
  list( APPEND physicsSolvers_sources
        python/PySolver.cpp )
  list( APPEND dependencyList Python3::Python pylvarray )
endif()

blt_add_library( NAME       physicsSolvers
                 SOURCES    ${physicsSolvers_sources}
                 HEADERS    ${physicsSolvers_headers}
                 DEPENDS_ON ${dependencyList} ${externalComponentDeps}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS} )

target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
if( externalComponentDeps )
  target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()


