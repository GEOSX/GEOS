/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */



/**
 * @file CompositionalMultiphaseFVM.cpp
 */

#include "CompositionalMultiphaseFVM.hpp"

#include "common/MpiWrapper.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "dataRepository/Group.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/IsothermalCompositionalMultiphaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalCompositionalMultiphaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/StabilizedCompositionalMultiphaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/DissipationCompositionalMultiphaseFVMKernels.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

CompositionalMultiphaseFVM::CompositionalMultiphaseFVM( const string & name,
                                                        Group * const parent )
  :
  CompositionalMultiphaseBase( name, parent )
{
  registerWrapper( viewKeyStruct::useDBCString(), &m_dbcParams.useDBC ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Enable Dissipation-based continuation flux" );

  registerWrapper( viewKeyStruct::omegaDBCString(), &m_dbcParams.omega ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor by which DBC flux is multiplied" );

  registerWrapper( viewKeyStruct::continuationDBCString(), &m_dbcParams.continuation ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for enabling continuation parameter" );

  registerWrapper( viewKeyStruct::miscibleDBCString(), &m_dbcParams.miscible ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for enabling DBC formulation with/without miscibility" );

  registerWrapper( viewKeyStruct::kappaminDBCString(), &m_dbcParams.kappamin ).
    setApplyDefaultValue( 1e-20 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor that controls how much dissipation is kept in the system when continuation is used" );

  registerWrapper( viewKeyStruct::contMultiplierDBCString(), &m_dbcParams.contMultiplier ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor by which continuation parameter is changed every newton when DBC is used" );

  registerWrapper( viewKeyStruct::scalingTypeString(), &m_scalingType ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( ScalingType::Global ).
    setDescription( "Solution scaling type."
                    "Valid options:\n* " + EnumStrings< ScalingType >::concat( "\n* " ) );
}

void CompositionalMultiphaseFVM::postInputInitialization()
{
  CompositionalMultiphaseBase::postInputInitialization();

  if( m_scalingType == ScalingType::Local && m_nonlinearSolverParameters.m_lineSearchAction != NonlinearSolverParameters::LineSearchAction::None )
  {
    GEOS_ERROR( GEOS_FMT( "{}: line search is not supported for {} = {}", getName(), viewKeyStruct::scalingTypeString(), EnumStrings< ScalingType >::toString( ScalingType::Local )) );
  }
}

void CompositionalMultiphaseFVM::initializePreSubGroups()
{
  CompositionalMultiphaseBase::initializePreSubGroups();

  m_linearSolverParameters.get().mgr.strategy = m_isThermal
                                                ? LinearSolverParameters::MGR::StrategyType::thermalCompositionalMultiphaseFVM
                                                : LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseFVM;

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  if( !fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    GEOS_ERROR( "A discretization deriving from FluxApproximationBase must be selected with CompositionalMultiphaseFlow" );
  }

}

void CompositionalMultiphaseFVM::setupDofs( DomainPartition const & domain,
                                            DofManager & dofManager ) const
{
  // add a field for the cell-centered degrees of freedom
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       getMeshTargets() );

  // this call with instruct GEOS to reorder the dof numbers
  dofManager.setLocalReorderingType( viewKeyStruct::elemDofFieldString(),
                                     DofManager::LocalReorderingType::ReverseCutHillMcKee );

  // for the volume balance equation, disable global coupling
  // this equation is purely local (not coupled to neighbors or other physics)
  dofManager.disableGlobalCouplingForEquation( viewKeyStruct::elemDofFieldString(),
                                               m_numComponents );


  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), fluxApprox );
}


void CompositionalMultiphaseFVM::assembleFluxTerms( real64 const dt,
                                                    DomainPartition const & domain,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      // Convective flux
      if( m_isThermal )
      {
        thermalCompositionalMultiphaseFVMKernels::
          FaceBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     elemDofKey,
                                                     m_hasCapPressure,
                                                     m_useTotalMassEquation,
                                                     getName(),
                                                     mesh.getElemManager(),
                                                     stencilWrapper,
                                                     dt,
                                                     localMatrix.toViewConstSizes(),
                                                     localRhs.toView() );
      }
      else
      {
        if( m_dbcParams.useDBC )
        {
          dissipationCompositionalMultiphaseFVMKernels::
            FaceBasedAssemblyKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       m_numPhases,
                                                       dofManager.rankOffset(),
                                                       elemDofKey,
                                                       m_hasCapPressure,
                                                       m_useTotalMassEquation,
                                                       getName(),
                                                       mesh.getElemManager(),
                                                       stencilWrapper,
                                                       dt,
                                                       localMatrix.toViewConstSizes(),
                                                       localRhs.toView(),
                                                       m_dbcParams.omega,
                                                       getNonlinearSolverParameters().m_numNewtonIterations,
                                                       m_dbcParams.continuation,
                                                       m_dbcParams.miscible,
                                                       m_dbcParams.kappamin,
                                                       m_dbcParams.contMultiplier );
        }
        else
        {
          isothermalCompositionalMultiphaseFVMKernels::
            FaceBasedAssemblyKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       m_numPhases,
                                                       dofManager.rankOffset(),
                                                       elemDofKey,
                                                       m_hasCapPressure,
                                                       m_useTotalMassEquation,
                                                       fluxApprox.upwindingParams(),
                                                       getName(),
                                                       mesh.getElemManager(),
                                                       stencilWrapper,
                                                       dt,
                                                       localMatrix.toViewConstSizes(),
                                                       localRhs.toView() );
        }
      }

      // Diffusive and dispersive flux
      if( m_hasDiffusion || m_hasDispersion )
      {

        if( m_isThermal )
        {
          thermalCompositionalMultiphaseFVMKernels::
            DiffusionDispersionFaceBasedAssemblyKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       m_numPhases,
                                                       dofManager.rankOffset(),
                                                       elemDofKey,
                                                       m_hasDiffusion,
                                                       m_hasDispersion,
                                                       m_useTotalMassEquation,
                                                       getName(),
                                                       mesh.getElemManager(),
                                                       stencilWrapper,
                                                       dt,
                                                       localMatrix.toViewConstSizes(),
                                                       localRhs.toView() );
        }
        else
        {
          isothermalCompositionalMultiphaseFVMKernels::
            DiffusionDispersionFaceBasedAssemblyKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       m_numPhases,
                                                       dofManager.rankOffset(),
                                                       elemDofKey,
                                                       m_hasDiffusion,
                                                       m_hasDispersion,
                                                       m_useTotalMassEquation,
                                                       getName(),
                                                       mesh.getElemManager(),
                                                       stencilWrapper,
                                                       dt,
                                                       localMatrix.toViewConstSizes(),
                                                       localRhs.toView() );
        }
      }

    } );
  } );
}

void CompositionalMultiphaseFVM::assembleStabilizedFluxTerms( real64 const dt,
                                                              DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      // Thermal implementation not supported yet

      stabilizedCompositionalMultiphaseFVMKernels::
        FaceBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                   m_numPhases,
                                                   dofManager.rankOffset(),
                                                   elemDofKey,
                                                   m_hasCapPressure,
                                                   m_useTotalMassEquation,
                                                   getName(),
                                                   mesh.getElemManager(),
                                                   stencilWrapper,
                                                   dt,
                                                   localMatrix.toViewConstSizes(),
                                                   localRhs.toView() );

    } );
  } );
}

real64 CompositionalMultiphaseFVM::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                          real64 const & GEOS_UNUSED_PARAM( dt ),
                                                          DomainPartition const & domain,
                                                          DofManager const & dofManager,
                                                          arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  integer constexpr numNorm = 3; // mass/volume balance and energy balance
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );

  physicsSolverBaseKernels::NormType const normType = getNonlinearSolverParameters().normType();

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      real64 subRegionResidualNorm[numNorm]{};
      real64 subRegionResidualNormalizer[numNorm]{};

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      // step 1: compute the norm in the subRegion

      if( m_isThermal )
      {
        string const & solidInternalEnergyName = subRegion.getReference< string >( viewKeyStruct::solidInternalEnergyNamesString() );
        SolidInternalEnergy const & solidInternalEnergy = getConstitutiveModel< SolidInternalEnergy >( subRegion, solidInternalEnergyName );

        thermalCompositionalMultiphaseBaseKernels::
          ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( normType,
                                                     numFluidComponents(),
                                                     numFluidPhases(),
                                                     rankOffset,
                                                     dofKey,
                                                     localRhs,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     solidInternalEnergy,
                                                     m_nonlinearSolverParameters.m_minNormalizer,
                                                     subRegionResidualNorm,
                                                     subRegionResidualNormalizer );
      }
      else
      {
        real64 subRegionFlowResidualNorm[2]{};
        real64 subRegionFlowResidualNormalizer[2]{};
        isothermalCompositionalMultiphaseBaseKernels::
          ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( normType,
                                                     numFluidComponents(),
                                                     rankOffset,
                                                     dofKey,
                                                     localRhs,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     m_nonlinearSolverParameters.m_minNormalizer,
                                                     subRegionFlowResidualNorm,
                                                     subRegionFlowResidualNormalizer );
        // mass
        subRegionResidualNorm[0] = subRegionFlowResidualNorm[0];
        subRegionResidualNormalizer[0] = subRegionFlowResidualNormalizer[0];
        // volume
        subRegionResidualNorm[1] = subRegionFlowResidualNorm[1];
        subRegionResidualNormalizer[1] = subRegionFlowResidualNormalizer[1];
      }

      // step 2: first reduction across meshBodies/regions/subRegions

      if( normType == physicsSolverBaseKernels::NormType::Linf )
      {
        physicsSolverBaseKernels::LinfResidualNormHelper::
          updateLocalNorm< numNorm >( subRegionResidualNorm, localResidualNorm );
      }
      else
      {
        physicsSolverBaseKernels::L2ResidualNormHelper::
          updateLocalNorm< numNorm >( subRegionResidualNorm, subRegionResidualNormalizer, localResidualNorm, localResidualNormalizer );
      }
    } );
  } );

  // step 3: second reduction across MPI ranks

  real64 residualNorm = 0.0;
  if( m_isThermal )
  {
    array1d< real64 > globalResidualNorm;
    globalResidualNorm.resize( numNorm );
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      physicsSolverBaseKernels::LinfResidualNormHelper::
        computeGlobalNorm( localResidualNorm, globalResidualNorm );
    }
    else
    {
      physicsSolverBaseKernels::L2ResidualNormHelper::
        computeGlobalNorm( localResidualNorm, localResidualNormalizer, globalResidualNorm );
    }
    residualNorm = sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1]  + globalResidualNorm[2] * globalResidualNorm[2] );

    GEOS_LOG_LEVEL_INFO_RANK_0_NLR( logInfo::Convergence, GEOS_FMT( "        ( Rmass Rvol ) = ( {:4.2e} {:4.2e} )        ( Renergy ) = ( {:4.2e} )",
                                                                    globalResidualNorm[0], globalResidualNorm[1], globalResidualNorm[2] ));
  }
  else
  {
    array1d< real64 > globalResidualNorm;
    globalResidualNorm.resize( numNorm - 1 );
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      physicsSolverBaseKernels::LinfResidualNormHelper::
        computeGlobalNorm( localResidualNorm, globalResidualNorm );
    }
    else
    {
      physicsSolverBaseKernels::L2ResidualNormHelper::
        computeGlobalNorm( localResidualNorm, localResidualNormalizer, globalResidualNorm );
    }
    residualNorm = sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1] );

    if( getLogLevel() >= 1 && logger::internal::rank == 0 )
    {
      std::cout << GEOS_FMT( "        ( Rmass Rvol ) = ( {:4.2e} {:4.2e} )",
                             globalResidualNorm[0], globalResidualNorm[1] );
    }
  }

  return residualNorm;
}

real64 CompositionalMultiphaseFVM::scalingForSystemSolution( DomainPartition & domain,
                                                             DofManager const & dofManager,
                                                             arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  real64 scalingFactor = 1.0;
  real64 maxDeltaPres = 0.0, maxDeltaCompDens = 0.0, maxDeltaTemp = 0.0;
  real64 minPresScalingFactor = 1.0, minCompDensScalingFactor = 1.0, minTempScalingFactor = 1.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const pressure = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const temperature = subRegion.getField< fields::flow::temperature >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< fields::flow::globalCompDensity >();
      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
      arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< fields::flow::temperatureScalingFactor >();
      arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::flow::globalCompDensityScalingFactor >();

      const integer temperatureOffset = m_numComponents+1;
      auto const subRegionData =
        m_isThermal
  ? thermalCompositionalMultiphaseBaseKernels::
          ScalingForSystemSolutionKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                     m_maxAbsolutePresChange,
                                                     m_maxRelativeTempChange,
                                                     m_maxCompFracChange,
                                                     m_maxRelativeCompDensChange,
                                                     pressure,
                                                     temperature,
                                                     compDens,
                                                     pressureScalingFactor,
                                                     compDensScalingFactor,
                                                     temperatureScalingFactor,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     dofKey,
                                                     subRegion,
                                                     localSolution,
                                                     temperatureOffset )
  : isothermalCompositionalMultiphaseBaseKernels::
          ScalingForSystemSolutionKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                     m_maxAbsolutePresChange,
                                                     m_maxCompFracChange,
                                                     m_maxRelativeCompDensChange,
                                                     pressure,
                                                     compDens,
                                                     pressureScalingFactor,
                                                     compDensScalingFactor,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     dofKey,
                                                     subRegion,
                                                     localSolution );

      if( m_scalingType == ScalingType::Global )
      {
        scalingFactor = std::min( scalingFactor, subRegionData.localMinVal );
      }
      maxDeltaPres  = std::max( maxDeltaPres, subRegionData.localMaxDeltaPres );
      maxDeltaCompDens = std::max( maxDeltaCompDens, subRegionData.localMaxDeltaCompDens );
      maxDeltaTemp = std::max( maxDeltaTemp, subRegionData.localMaxDeltaTemp );
      minPresScalingFactor = std::min( minPresScalingFactor, subRegionData.localMinPresScalingFactor );
      minCompDensScalingFactor = std::min( minCompDensScalingFactor, subRegionData.localMinCompDensScalingFactor );
      minTempScalingFactor = std::min( minTempScalingFactor, subRegionData.localMinTempScalingFactor );
    } );
  } );

  scalingFactor = MpiWrapper::min( scalingFactor );
  maxDeltaPres  = MpiWrapper::max( maxDeltaPres );
  maxDeltaCompDens = MpiWrapper::max( maxDeltaCompDens );
  minPresScalingFactor = MpiWrapper::min( minPresScalingFactor );
  minCompDensScalingFactor = MpiWrapper::min( minCompDensScalingFactor );

  string const massUnit = m_useMass ? "kg/m3" : "mol/m3";
  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Max pressure change = {} Pa (before scaling)",
                                                           getName(), GEOS_FMT( "{:.{}f}", maxDeltaPres, 3 ) ) );
  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Max component density change = {} {} (before scaling)",
                                                           getName(), GEOS_FMT( "{:.{}f}", maxDeltaCompDens, 3 ), massUnit ) );

  if( m_isThermal )
  {
    maxDeltaTemp = MpiWrapper::max( maxDeltaTemp );
    minTempScalingFactor = MpiWrapper::min( minTempScalingFactor );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Max temperature change = {} K (before scaling)",
                                                             getName(), GEOS_FMT( "{:.{}f}", maxDeltaTemp, 3 ) ) );
  }

  if( m_scalingType == ScalingType::Local )
  {
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Min pressure scaling factor = {}", getName(), minPresScalingFactor ) );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Min component density scaling factor = {}", getName(), minCompDensScalingFactor ) );
    if( m_isThermal )
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Min temperature scaling factor = {}", getName(), minTempScalingFactor ) );
    }
  }

  return LvArray::math::max( scalingFactor, m_minScalingFactor );
}

bool CompositionalMultiphaseFVM::checkSystemSolution( DomainPartition & domain,
                                                      DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localSolution,
                                                      real64 const scalingFactor )
{
  GEOS_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  integer localCheck = 1;
  real64 minPres = 0.0, minDens = 0.0, minTotalDens = 0.0;
  integer numNegPres = 0, numNegDens = 0, numNegTotalDens = 0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const pressure =
        subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const temperature =
        subRegion.getField< fields::flow::temperature >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens =
        subRegion.getField< fields::flow::globalCompDensity >();
      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
      arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< fields::flow::temperatureScalingFactor >();
      arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::flow::globalCompDensityScalingFactor >();
      // check that pressure and component densities are non-negative
      // for thermal, check that temperature is above 273.15 K
      const integer temperatureOffset = m_numComponents+1;
      auto const subRegionData =
        m_isThermal
  ? thermalCompositionalMultiphaseBaseKernels::
          SolutionCheckKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_allowCompDensChopping,
                                                     m_allowNegativePressure,
                                                     m_scalingType,
                                                     scalingFactor,
                                                     pressure,
                                                     temperature,
                                                     compDens,
                                                     pressureScalingFactor,
                                                     temperatureScalingFactor,
                                                     compDensScalingFactor,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     dofKey,
                                                     subRegion,
                                                     localSolution,
                                                     temperatureOffset )
  : isothermalCompositionalMultiphaseBaseKernels::
          SolutionCheckKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_allowCompDensChopping,
                                                     m_allowNegativePressure,
                                                     m_scalingType,
                                                     scalingFactor,
                                                     pressure,

                                                     compDens,
                                                     pressureScalingFactor,

                                                     compDensScalingFactor,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     dofKey,
                                                     subRegion,
                                                     localSolution );

      localCheck = std::min( localCheck, subRegionData.localMinVal );

      minPres  = std::min( minPres, subRegionData.localMinPres );
      minDens = std::min( minDens, subRegionData.localMinDens );
      minTotalDens = std::min( minTotalDens, subRegionData.localMinTotalDens );
      numNegPres += subRegionData.localNumNegPressures;
      numNegDens += subRegionData.localNumNegDens;
      numNegTotalDens += subRegionData.localNumNegTotalDens;
    } );
  } );

  minPres  = MpiWrapper::min( minPres );
  minDens = MpiWrapper::min( minDens );
  minTotalDens = MpiWrapper::min( minTotalDens );
  numNegPres = MpiWrapper::sum( numNegPres );
  numNegDens = MpiWrapper::sum( numNegDens );
  numNegTotalDens = MpiWrapper::sum( numNegTotalDens );

  if( numNegPres > 0 )
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Number of negative pressure values: {}, minimum value: {} Pa",
                                                             getName(), numNegPres, fmt::format( "{:.{}f}", minPres, 3 ) ) );
  string const massUnit = m_useMass ? "kg/m3" : "mol/m3";
  if( numNegDens > 0 )
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Number of negative component density values: {}, minimum value: {} {}}",
                                                             getName(), numNegDens, fmt::format( "{:.{}f}", minDens, 3 ), massUnit ) );
  if( minTotalDens > 0 )
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Number of negative total density values: {}, minimum value: {} {}}",
                                                             getName(), minTotalDens, fmt::format( "{:.{}f}", minDens, 3 ), massUnit ) );

  return MpiWrapper::min( localCheck );
}

void CompositionalMultiphaseFVM::applySystemSolution( DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localSolution,
                                                      real64 const scalingFactor,
                                                      real64 const GEOS_UNUSED_PARAM( dt ),
                                                      DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  bool const localScaling = m_scalingType == ScalingType::Local;

  DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );
  DofManager::CompMask componentMask( m_numDofPerCell, 1, m_numComponents+1 );

  if( localScaling )
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 fields::flow::pressureScalingFactor::key(),
                                 pressureMask );
  }
  else
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 scalingFactor,
                                 pressureMask );
  }

  if( localScaling )
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::globalCompDensity::key(),
                                 fields::flow::globalCompDensityScalingFactor::key(),
                                 componentMask );
  }
  else
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::globalCompDensity::key(),
                                 scalingFactor,
                                 componentMask );
  }

  if( m_isThermal )
  {
    DofManager::CompMask temperatureMask( m_numDofPerCell, m_numComponents+1, m_numComponents+2 );
    if( localScaling )
    {
      dofManager.addVectorToField( localSolution,
                                   viewKeyStruct::elemDofFieldString(),
                                   fields::flow::temperature::key(),
                                   fields::flow::temperatureScalingFactor::key(),
                                   temperatureMask );
    }
    else
    {
      dofManager.addVectorToField( localSolution,
                                   viewKeyStruct::elemDofFieldString(),
                                   fields::flow::temperature::key(),
                                   scalingFactor,
                                   temperatureMask );
    }
  }

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    chopNegativeDensities( domain );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    std::vector< string > fields{ fields::flow::pressure::key(), fields::flow::globalCompDensity::key() };
    if( m_isThermal )
    {
      fields.emplace_back( fields::flow::temperature::key() );
    }
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( fields, regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}

void CompositionalMultiphaseFVM::updatePhaseMobility( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  // note that the phase mobility computed here also includes phase density
  string const & fluidName = dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, fluidName );

  string const & relpermName = dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() );
  RelativePermeabilityBase const & relperm = getConstitutiveModel< RelativePermeabilityBase >( dataGroup, relpermName );

  if( m_isThermal )
  {
    thermalCompositionalMultiphaseFVMKernels::
      PhaseMobilityKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 dataGroup,
                                                 fluid,
                                                 relperm );
  }
  else
  {
    isothermalCompositionalMultiphaseFVMKernels::
      PhaseMobilityKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 dataGroup,
                                                 fluid,
                                                 relperm );
  }
}

void CompositionalMultiphaseFVM::applyBoundaryConditions( real64 time_n,
                                                          real64 dt,
                                                          DomainPartition & domain,
                                                          DofManager const & dofManager,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  CompositionalMultiphaseBase::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
  if( !m_keepVariablesConstantDuringInitStep )
  {
    applyFaceDirichletBC( time_n, dt, dofManager, domain, localMatrix, localRhs );
  }
}

bool CompositionalMultiphaseFVM::validateFaceDirichletBC( DomainPartition & domain,
                                                          real64 const time ) const
{
  constexpr integer MAX_NC = MultiFluidBase::MAX_NUM_COMPONENTS;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  bool bcConsistent = true;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {

    // maps to check consistent application of BC
    // maps: setName (-> numComps)
    map< string, ComponentMask< MAX_NC > > bcPresCompStatusMap; // check that pressure/comp are present/consistent
    set< string > bcTempStatusMap; // check that temperature is present/consistent

    // 1. Check pressure Dirichlet BCs
    fsManager.apply< FaceManager >( time,
                                    mesh,
                                    fields::flow::pressure::key(),
                                    [&]( FieldSpecificationBase const &,
                                         string const & setName,
                                         SortedArrayView< localIndex const > const &,
                                         FaceManager &,
                                         string const & )
    {
      // Check whether pressure has already been applied to this set
      if( bcPresCompStatusMap.count( setName ) > 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Conflicting pressure boundary conditions on set {}", setName ) );
      }
      bcPresCompStatusMap[setName].setNumComp( m_numComponents );
    } );

    // 2. Check temperature Dirichlet BCs (we always require a temperature for face-based BCs)
    fsManager.apply< FaceManager >( time,
                                    mesh,
                                    fields::flow::temperature::key(),
                                    [&]( FieldSpecificationBase const &,
                                         string const & setName,
                                         SortedArrayView< localIndex const > const &,
                                         FaceManager &,
                                         string const & )
    {
      // 2.1 Check whether temperature has already been applied to this set
      if( bcTempStatusMap.count( setName ) > 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Conflicting temperature boundary conditions on set {}", setName ) );
      }
      bcTempStatusMap.insert( setName );

      // 2.2 Check that there is pressure bc applied to this set
      if( bcPresCompStatusMap.count( setName ) == 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Pressure boundary condition not prescribed on set {}", setName ) );
      }
    } );

    // 3. Check composition BC (global component fraction)
    fsManager.apply< FaceManager >( time,
                                    mesh,
                                    fields::flow::globalCompFraction::key(),
                                    [&] ( FieldSpecificationBase const & fs,
                                          string const & setName,
                                          SortedArrayView< localIndex const > const &,
                                          FaceManager &,
                                          string const & )
    {
      // 3.1 Check pressure, temperature, and record composition bc application
      integer const comp = fs.getComponent();

      if( bcPresCompStatusMap.count( setName ) == 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Pressure boundary condition not prescribed on set {}", setName ) );
      }
      if( bcTempStatusMap.count( setName ) == 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Temperature boundary condition not prescribed on set {}. \n"
                                "Note that for face boundary conditions, you must provide a temperature", setName ) );
      }
      if( comp < 0 || comp >= m_numComponents )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Invalid component index [{}] in composition boundary condition {}", comp, fs.getName() ) );
        return; // can't check next part with invalid component id
      }

      ComponentMask< MAX_NC > & compMask = bcPresCompStatusMap[setName];
      if( compMask[comp] )
      {
        bcConsistent = false;
        GEOS_WARNING( GEOS_FMT( "Conflicting composition[{}] boundary conditions on set {}", comp, setName ) );
      }
      compMask.set( comp );
    } );

    // 3.2 Check consistency between composition BC applied to sets
    for( auto const & setEntry : bcPresCompStatusMap )
    {
      ComponentMask< MAX_NC > const & compMask = setEntry.second;
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        if( !compMask[ic] )
        {
          bcConsistent = false;
          GEOS_WARNING( GEOS_FMT( "Boundary condition not applied to composition[{}] on set {}", ic, setEntry.first ) );
        }
      }
    }
  } );

  return bcConsistent;
}

namespace
{
char const faceBcLogMessage[] =
  "CompositionalMultiphaseFVM {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target faces (including ghost faces) is {}."
  "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set.";
}

void CompositionalMultiphaseFVM::applyFaceDirichletBC( real64 const time_n,
                                                       real64 const dt,
                                                       DofManager const & dofManager,
                                                       DomainPartition & domain,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // Only validate BC at the beginning of Newton loop
  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    bool const bcConsistent = validateFaceDirichletBC( domain, time_n + dt );
    GEOS_ERROR_IF( !bcConsistent, GEOS_FMT( "CompositionalMultiphaseBase {}: inconsistent boundary conditions", getDataContext() ) );
  }

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    // Take BCs defined for "pressure" field and apply values to "facePressure"
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    fields::flow::pressure::key(), fields::flow::facePressure::key() );
    // Take BCs defined for "globalCompFraction" field and apply values to "faceGlobalCompFraction"
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    fields::flow::globalCompFraction::key(), fields::flow::faceGlobalCompFraction::key() );
    // Take BCs defined for "temperature" field and apply values to "faceTemperature"
    applyFieldValue< FaceManager >( time_n, dt, mesh, faceBcLogMessage,
                                    fields::flow::temperature::key(), fields::flow::faceTemperature::key() );

    // Then launch the face Dirichlet kernel
    fsManager.apply< FaceManager >( time_n + dt,
                                    mesh,
                                    fields::flow::pressure::key(), // we have required that pressure is always present
                                    [&] ( FieldSpecificationBase const &,
                                          string const & setName,
                                          SortedArrayView< localIndex const > const &,
                                          FaceManager &,
                                          string const & )
    {
      BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
      if( stencil.size() == 0 )
      {
        return;
      }

      // TODO: same issue as in the single-phase case
      //       currently we just use model from the first cell in this stencil
      //       since it's not clear how to create fluid kernel wrappers for arbitrary models.
      //       Can we just use cell properties for an approximate flux computation?
      //       Then we can forget about capturing the fluid model.
      localIndex const er = stencil.getElementRegionIndices()( 0, 0 );
      localIndex const esr = stencil.getElementSubRegionIndices()( 0, 0 );
      ElementSubRegionBase & subRegion = elemManager.getRegion( er ).getSubRegion( esr );
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & multiFluidBase = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

      BoundaryStencilWrapper const stencilWrapper = stencil.createKernelWrapper();

      string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

      if( m_isThermal )
      {
        //todo (jafranc) extend upwindScheme name if satisfied in isothermalCase
        thermalCompositionalMultiphaseFVMKernels::
          DirichletFaceBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     m_useTotalMassEquation,
                                                     elemDofKey,
                                                     getName(),
                                                     faceManager,
                                                     elemManager,
                                                     stencilWrapper,
                                                     multiFluidBase,
                                                     dt,
                                                     localMatrix,
                                                     localRhs );
      }
      else
      {
        isothermalCompositionalMultiphaseFVMKernels::
          DirichletFaceBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_numPhases,
                                                     dofManager.rankOffset(),
                                                     m_useTotalMassEquation,
                                                     elemDofKey,
                                                     getName(),
                                                     faceManager,
                                                     elemManager,
                                                     stencilWrapper,
                                                     multiFluidBase,
                                                     dt,
                                                     localMatrix,
                                                     localRhs );
      }

    } );
  } );
}

void CompositionalMultiphaseFVM::applyAquiferBC( real64 const time,
                                                 real64 const dt,
                                                 DofManager const & dofManager,
                                                 DomainPartition & domain,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

    isothermalCompositionalMultiphaseFVMKernels::
      AquiferBCKernel::CompFlowAccessors compFlowAccessors( mesh.getElemManager(), getName() );
    isothermalCompositionalMultiphaseFVMKernels::
      AquiferBCKernel::MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), getName() );

    fsManager.apply< FaceManager,
                     AquiferBoundaryCondition >( time + dt,
                                                 mesh,
                                                 AquiferBoundaryCondition::catalogName(),
                                                 [&] ( AquiferBoundaryCondition const & bc,
                                                       string const & setName,
                                                       SortedArrayView< localIndex const > const &,
                                                       FaceManager & faceManager,
                                                       string const & )
    {
      BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
      if( bc.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( stencil.size() );
        GEOS_LOG_RANK_0( GEOS_FMT( faceBcLogMessage,
                                   getName(), time+dt, bc.getCatalogName(), bc.getName(),
                                   setName, faceManager.getName(), bc.getScale(), numTargetFaces ) );
      }

      if( stencil.size() == 0 )
      {
        return;
      }

      AquiferBoundaryCondition::KernelWrapper aquiferBCWrapper = bc.createKernelWrapper();
      bool const allowAllPhasesIntoAquifer = bc.allowAllPhasesIntoAquifer();
      localIndex const waterPhaseIndex = bc.getWaterPhaseIndex();
      real64 const & aquiferWaterPhaseDens = bc.getWaterPhaseDensity();
      arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac = bc.getWaterPhaseComponentFraction();

      // While this kernel is waiting for a factory class, pass all the accessors here
      isothermalCompositionalMultiphaseBaseKernels::KernelLaunchSelector1
      < isothermalCompositionalMultiphaseFVMKernels::AquiferBCKernel >( m_numComponents,
                                                                        m_numPhases,
                                                                        waterPhaseIndex,
                                                                        allowAllPhasesIntoAquifer,
                                                                        m_useTotalMassEquation,
                                                                        stencil,
                                                                        dofManager.rankOffset(),
                                                                        elemDofNumber.toNestedViewConst(),
                                                                        aquiferBCWrapper,
                                                                        aquiferWaterPhaseDens,
                                                                        aquiferWaterPhaseCompFrac,
                                                                        compFlowAccessors.get( fields::ghostRank{} ),
                                                                        compFlowAccessors.get( fields::flow::pressure{} ),
                                                                        compFlowAccessors.get( fields::flow::pressure_n{} ),
                                                                        compFlowAccessors.get( fields::flow::gravityCoefficient{} ),
                                                                        compFlowAccessors.get( fields::flow::phaseVolumeFraction{} ),
                                                                        compFlowAccessors.get( fields::flow::dPhaseVolumeFraction{} ),
                                                                        compFlowAccessors.get( fields::flow::dGlobalCompFraction_dGlobalCompDensity{} ),
                                                                        multiFluidAccessors.get( fields::multifluid::phaseDensity{} ),
                                                                        multiFluidAccessors.get( fields::multifluid::dPhaseDensity{} ),
                                                                        multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} ),
                                                                        multiFluidAccessors.get( fields::multifluid::dPhaseCompFraction{} ),
                                                                        time,
                                                                        dt,
                                                                        localMatrix.toViewConstSizes(),
                                                                        localRhs.toView() );
    } );
  } );

}

//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, CompositionalMultiphaseFVM, string const &, Group * const )
//END_SPHINX_INCLUDE_01
}// namespace geos
