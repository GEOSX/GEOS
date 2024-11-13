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
 * @file SinglePhaseFVM.cpp
 */

#include "SinglePhaseFVM.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidSelector.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/AccumulationKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/ThermalAccumulationKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/ResidualNormKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/ThermalFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/DirichletFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/ThermalDirichletFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/StabilizedFluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/AquiferBCKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/proppant/ProppantFluxKernels.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsEmbeddedFractures.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanicsEmbeddedFractures.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsConformingFractures.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanicsConformingFractures.hpp"

/**
 * @namespace the geos namespace that encapsulates the majority of the code
 */
namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace singlePhaseBaseKernels;
using namespace singlePhaseFVMKernels;

template< typename BASE >
SinglePhaseFVM< BASE >::SinglePhaseFVM( const string & name,
                                        Group * const parent ):
  BASE( name, parent )
{}

template< typename BASE >
void SinglePhaseFVM< BASE >::initializePreSubGroups()
{
  BASE::initializePreSubGroups();

  DomainPartition & domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();

  if( !fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    GEOS_ERROR( "A discretization deriving from FluxApproximationBase must be selected with SinglePhaseFVM" );
  }

  if( m_isThermal )
  {
    // For thermal simulations 2 pdes are considered so we let AMG know.
    LinearSolverParameters & linParams = m_linearSolverParameters.get();
    linParams.amg.numFunctions = 2;
  }
}

template< typename BASE >
void SinglePhaseFVM< BASE >::setupDofs( DomainPartition const & domain,
                                        DofManager & dofManager ) const
{
  dofManager.addField( BASE::viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       BASE::getMeshTargets() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( BASE::viewKeyStruct::elemDofFieldString(), fluxApprox );
}

template< typename BASE >
void SinglePhaseFVM< BASE >::setupSystem( DomainPartition & domain,
                                          DofManager & dofManager,
                                          CRSMatrix< real64, globalIndex > & localMatrix,
                                          ParallelVector & rhs,
                                          ParallelVector & solution,
                                          bool const setSparsity )
{
  GEOS_MARK_FUNCTION;
  BASE::setupSystem( domain,
                     dofManager,
                     localMatrix,
                     rhs,
                     solution,
                     setSparsity );

}

template< typename BASE >
real64 SinglePhaseFVM< BASE >::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                      real64 const & GEOS_UNUSED_PARAM( dt ),
                                                      DomainPartition const & domain,
                                                      DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  integer constexpr numNorm = 2; // mass balance and energy balance
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );

  physicsSolverBaseKernels::NormType const normType = BASE::getNonlinearSolverParameters().normType();

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( BASE::viewKeyStruct::elemDofFieldString() );

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      real64 subRegionResidualNorm[numNorm]{};
      real64 subRegionResidualNormalizer[numNorm]{};

      // step 1: compute the norm in the subRegion

      if( m_isThermal )
      {
        singlePhaseBaseKernels::
          ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( normType,
                                                     rankOffset,
                                                     dofKey,
                                                     localRhs,
                                                     subRegion,
                                                     m_nonlinearSolverParameters.m_minNormalizer,
                                                     subRegionResidualNorm,
                                                     subRegionResidualNormalizer );
      }
      else
      {
        real64 subRegionFlowResidualNorm[1]{};
        real64 subRegionFlowResidualNormalizer[1]{};
        singlePhaseBaseKernels::
          ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( normType,
                                                     rankOffset,
                                                     dofKey,
                                                     localRhs,
                                                     subRegion,
                                                     m_nonlinearSolverParameters.m_minNormalizer,
                                                     subRegionFlowResidualNorm,
                                                     subRegionFlowResidualNormalizer );
        subRegionResidualNorm[0] = subRegionFlowResidualNorm[0];
        subRegionResidualNormalizer[0] = subRegionFlowResidualNormalizer[0];
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
    residualNorm = sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1] );

    GEOS_LOG_LEVEL_INFO_RANK_0_NLR( logInfo::Convergence, GEOS_FMT( "        ( R{} ) = ( {:4.2e} )        ( Renergy ) = ( {:4.2e} )",
                                                                    FlowSolverBase::coupledSolverAttributePrefix(), globalResidualNorm[0], globalResidualNorm[1] ));
  }
  else
  {

    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      physicsSolverBaseKernels::LinfResidualNormHelper::computeGlobalNorm( localResidualNorm[0], residualNorm );
    }
    else
    {
      physicsSolverBaseKernels::L2ResidualNormHelper::computeGlobalNorm( localResidualNorm[0], localResidualNormalizer[0], residualNorm );
    }

    GEOS_LOG_LEVEL_INFO_RANK_0_NLR( logInfo::Convergence,
                                    GEOS_FMT( "        ( R{} ) = ( {:4.2e} )", FlowSolverBase::coupledSolverAttributePrefix(), residualNorm ));
  }
  return residualNorm;
}


template< typename BASE >
void SinglePhaseFVM< BASE >::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  real64 const dt,
                                                  DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );
  if( m_isThermal )
  {
    DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );
    DofManager::CompMask temperatureMask( m_numDofPerCell, 1, 2 );

    dofManager.addVectorToField( localSolution,
                                 BASE::viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 scalingFactor,
                                 pressureMask );

    dofManager.addVectorToField( localSolution,
                                 BASE::viewKeyStruct::elemDofFieldString(),
                                 fields::flow::temperature::key(),
                                 scalingFactor,
                                 temperatureMask );
  }
  else
  {
    dofManager.addVectorToField( localSolution,
                                 BASE::viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 scalingFactor );
  }

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      arrayView1d< string const > const & regionNames )
  {
    std::vector< string > fields{ fields::flow::pressure::key() };

    if( m_isThermal )
    {
      fields.emplace_back( fields::flow::temperature::key() );
    }

    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( fields, regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors() );
  } );
}

template<>
void SinglePhaseFVM<>::assembleFluxTerms( real64 const dt,
                                          DomainPartition const & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & )
  {
    fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();


      if( m_isThermal )
      {
        thermalSinglePhaseFVMKernels::
          FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                               dofKey,
                                                                               getName(),
                                                                               mesh.getElemManager(),
                                                                               stencilWrapper,
                                                                               dt,
                                                                               localMatrix.toViewConstSizes(),
                                                                               localRhs.toView() );
      }
      else
      {
        singlePhaseFVMKernels::
          FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                               dofKey,
                                                                               getName(),
                                                                               mesh.getElemManager(),
                                                                               stencilWrapper,
                                                                               dt,
                                                                               localMatrix.toViewConstSizes(),
                                                                               localRhs.toView() );
      }


    } );
  } );

}


template< >
void SinglePhaseFVM< SinglePhaseBase >::assembleStabilizedFluxTerms( real64 const dt,
                                                                     DomainPartition const & domain,
                                                                     DofManager const & dofManager,
                                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                     arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & )
  {
    fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      // No thermal support yet
      stabilizedSinglePhaseFVMKernels::FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                                                            dofKey,
                                                                                                            getName(),
                                                                                                            mesh.getElemManager(),
                                                                                                            stencilWrapper,
                                                                                                            dt,
                                                                                                            localMatrix.toViewConstSizes(),
                                                                                                            localRhs.toView() );

    } );
  } );

}



template<>
void SinglePhaseFVM< SinglePhaseProppantBase >::assembleFluxTerms( real64 const dt,
                                                                   DomainPartition const & domain,
                                                                   DofManager const & dofManager,
                                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                   arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
    elemDofNumber = elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      typename FluxComputeKernelBase::SinglePhaseFlowAccessors flowAccessors( elemManager, getName() );
      typename FluxComputeKernelBase::SlurryFluidAccessors fluidAccessors( elemManager, getName() );
      typename FluxComputeKernelBase::ProppantPermeabilityAccessors permAccessors( elemManager, getName() );

      singlePhaseProppantFluxKernels::FaceElementFluxKernel::launch( stencilWrapper,
                                                                     dt,
                                                                     dofManager.rankOffset(),
                                                                     elemDofNumber.toNestedViewConst(),
                                                                     flowAccessors.get< fields::ghostRank >(),
                                                                     flowAccessors.get< fields::flow::pressure >(),
                                                                     flowAccessors.get< fields::flow::gravityCoefficient >(),
                                                                     fluidAccessors.get< fields::singlefluid::density >(),
                                                                     fluidAccessors.get< fields::singlefluid::dDensity_dPressure >(),
                                                                     flowAccessors.get< fields::flow::mobility >(),
                                                                     flowAccessors.get< fields::flow::dMobility_dPressure >(),
                                                                     permAccessors.get< fields::permeability::permeability >(),
                                                                     permAccessors.get< fields::permeability::dPerm_dPressure >(),
                                                                     permAccessors.get< fields::permeability::dPerm_dDispJump >(),
                                                                     permAccessors.get< fields::permeability::permeabilityMultiplier >(),
                                                                     this->gravityVector(),
                                                                     localMatrix,
                                                                     localRhs );
    } );
  } );
}

template< >
void SinglePhaseFVM< SinglePhaseProppantBase >::assembleStabilizedFluxTerms( real64 const dt,
                                                                             DomainPartition const & domain,
                                                                             DofManager const & dofManager,
                                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                             arrayView1d< real64 > const & localRhs )
{
  GEOS_UNUSED_VAR( dt, domain, dofManager, localMatrix, localRhs );
  GEOS_ERROR( "Stabilized flux not available with this flow solver" );
}

template< typename BASE >
void SinglePhaseFVM< BASE >::assembleEDFMFluxTerms( real64 const GEOS_UNUSED_PARAM ( time_n ),
                                                    real64 const dt,
                                                    DomainPartition const & domain,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs,
                                                    string const & jumpDofKey )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & )
  {
    fluxApprox.forStencils< CellElementStencilTPFA, EmbeddedSurfaceToCellStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      if( m_isThermal )
      {
        thermalSinglePhaseFVMKernels::
          FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                               dofKey,
                                                                               this->getName(),
                                                                               mesh.getElemManager(),
                                                                               stencilWrapper,
                                                                               dt,
                                                                               localMatrix.toViewConstSizes(),
                                                                               localRhs.toView() );
      }
      else
      {
        singlePhaseFVMKernels::
          FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                               dofKey,
                                                                               this->getName(),
                                                                               mesh.getElemManager(),
                                                                               stencilWrapper,
                                                                               dt,
                                                                               localMatrix.toViewConstSizes(),
                                                                               localRhs.toView() );
      }
    } );

    // Eventually, EmbeddedSurfaceToCellStencil should be moved here to account for the permeability of the fracture in the matrix-frac
    // connection
    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      if( m_isThermal )
      {
        thermalSinglePhasePoromechanicsEmbeddedFracturesKernels::
          ConnectorBasedAssemblyKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                                          dofKey,
                                                                                          jumpDofKey,
                                                                                          this->getName(),
                                                                                          mesh.getElemManager(),
                                                                                          stencilWrapper,
                                                                                          dt,
                                                                                          localMatrix.toViewConstSizes(),
                                                                                          localRhs.toView() );
      }
      else
      {
        singlePhasePoromechanicsEmbeddedFracturesKernels::
          ConnectorBasedAssemblyKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                                          dofKey,
                                                                                          jumpDofKey,
                                                                                          this->getName(),
                                                                                          mesh.getElemManager(),
                                                                                          stencilWrapper,
                                                                                          dt,
                                                                                          localMatrix.toViewConstSizes(),
                                                                                          localRhs.toView() );
      }
    } );
  } );

}

template< typename BASE >
void SinglePhaseFVM< BASE >::assembleHydrofracFluxTerms( real64 const GEOS_UNUSED_PARAM ( time_n ),
                                                         real64 const dt,
                                                         DomainPartition const & domain,
                                                         DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs,
                                                         CRSMatrixView< real64, localIndex const > const & dR_dAper )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );


  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & )
  {
    fluxApprox.forStencils< CellElementStencilTPFA, FaceElementToCellStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      if( m_isThermal )
      {
        thermalSinglePhaseFVMKernels::
          FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                               dofKey,
                                                                               this->getName(),
                                                                               mesh.getElemManager(),
                                                                               stencilWrapper,
                                                                               dt,
                                                                               localMatrix.toViewConstSizes(),
                                                                               localRhs.toView() );
      }
      else
      {
        singlePhaseFVMKernels::
          FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                               dofKey,
                                                                               this->getName(),
                                                                               mesh.getElemManager(),
                                                                               stencilWrapper,
                                                                               dt,
                                                                               localMatrix.toViewConstSizes(),
                                                                               localRhs.toView() );
      }
    } );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      if( m_isThermal )
      {
        thermalSinglePhasePoromechanicsConformingFracturesKernels::
          ConnectorBasedAssemblyKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                                          dofKey,
                                                                                          this->getName(),
                                                                                          mesh.getElemManager(),
                                                                                          stencilWrapper,
                                                                                          dt,
                                                                                          localMatrix.toViewConstSizes(),
                                                                                          localRhs.toView(),
                                                                                          dR_dAper );
      }
      else
      {
        singlePhasePoromechanicsConformingFracturesKernels::
          ConnectorBasedAssemblyKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                                          dofKey,
                                                                                          this->getName(),
                                                                                          mesh.getElemManager(),
                                                                                          stencilWrapper,
                                                                                          dt,
                                                                                          localMatrix.toViewConstSizes(),
                                                                                          localRhs.toView(),
                                                                                          dR_dAper );
      }
    } );
  } );


}

template< typename BASE >
void
SinglePhaseFVM< BASE >::applyBoundaryConditions( real64 const time_n,
                                                 real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  BASE::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
  if( !BASE::m_keepVariablesConstantDuringInitStep )
  {
    applyFaceDirichletBC( time_n, dt, dofManager, domain, localMatrix, localRhs );
  }
}

namespace
{
char const faceBcLogMessage[] =
  "SinglePhaseFVM {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. "
  "\nThe total number of target faces (including ghost faces) is {}. "
  "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set.";

GEOS_MAYBE_UNUSED char const incompleteBCLogmessage[] = "SinglePhaseFVM {}: at time {}, one or more Face boundary conditions are not complete. "
                                                        "Both pressure and temperature must be specified, one is missing.";

}


template< typename BASE >
void SinglePhaseFVM< BASE >::applyFaceDirichletBC( real64 const time_n,
                                                   real64 const dt,
                                                   DofManager const & dofManager,
                                                   DomainPartition & domain,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( BASE::viewKeyStruct::elemDofFieldString() );

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      arrayView1d< string const > const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    if( m_isThermal )
    {
      std::set< string > pressureSets;
      std::set< string > temperatureSets;
      // Take BCs defined for "pressure" field and apply values to "facePressure"
      fsManager.apply< FaceManager >( time_n + dt,
                                      mesh,
                                      fields::flow::pressure::key(),
                                      [&] ( FieldSpecificationBase const & fs,
                                            string const & setName,
                                            SortedArrayView< localIndex const > const & targetSet,
                                            FaceManager & targetGroup,
                                            string const & )
      {
        BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );

        if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
        {
          globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( stencil.size() );
          GEOS_LOG_RANK_0( GEOS_FMT( faceBcLogMessage,
                                     this->getName(), time_n+dt, fs.getCatalogName(), fs.getName(),
                                     setName, targetGroup.getName(), numTargetFaces ) );
        }

        if( stencil.size() == 0 )
        {
          return;
        }

        pressureSets.insert( setName );
        // first, evaluate BC to get primary field values (pressure)
        fs.applyFieldValue< FieldSpecificationEqual,
                            parallelDevicePolicy<> >( targetSet,
                                                      time_n + dt,
                                                      targetGroup,
                                                      fields::flow::facePressure::key() );
      } );

      // Take BCs defined for "temperature" field and apply values to "faceTemperature"
      fsManager.apply< FaceManager >( time_n + dt,
                                      mesh,
                                      fields::flow::temperature::key(),
                                      [&] ( FieldSpecificationBase const & fs,
                                            string const & setName,
                                            SortedArrayView< localIndex const > const & targetSet,
                                            FaceManager & targetGroup,
                                            string const & )
      {
        BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );

        if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
        {
          globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( stencil.size() );
          GEOS_LOG_RANK_0( GEOS_FMT( faceBcLogMessage,
                                     this->getName(), time_n+dt, fs.getCatalogName(), fs.getName(),
                                     setName, targetGroup.getName(), numTargetFaces ) );
        }

        if( stencil.size() == 0 )
        {
          return;
        }

        temperatureSets.insert( setName );
        // Specify the bc value of the field
        fs.applyFieldValue< FieldSpecificationEqual,
                            parallelDevicePolicy<> >( targetSet,
                                                      time_n + dt,
                                                      targetGroup,
                                                      fields::flow::faceTemperature::key() );

      } );

      GEOS_ERROR_IF( pressureSets != temperatureSets, GEOS_FMT( incompleteBCLogmessage, this->getName(), time_n + dt ) );

      // Take BCs defined for "temperature" field and apply values to "faceTemperature"
      for( auto const & setName : temperatureSets )
      {
        BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
        BoundaryStencilWrapper const stencilWrapper = stencil.createKernelWrapper();

        // TODO: currently we just use model from the first cell in this stencil
        //       since it's not clear how to create fluid kernel wrappers for arbitrary models.
        //       Can we just use cell properties for an approximate flux computation?
        //       Then we can forget about capturing the fluid model.
        localIndex const er = stencil.getElementRegionIndices()( 0, 0 );
        localIndex const esr = stencil.getElementSubRegionIndices()( 0, 0 );
        ElementSubRegionBase & subRegion = mesh.getElemManager().getRegion( er ).getSubRegion( esr );
        string const & fluidName = subRegion.getReference< string >( BASE::viewKeyStruct::fluidNamesString() );
        SingleFluidBase & fluidBase = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );

        thermalSinglePhaseFVMKernels::
          DirichletFluxComputeKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                     dofKey,
                                                     this->getName(),
                                                     faceManager,
                                                     elemManager,
                                                     stencilWrapper,
                                                     fluidBase,
                                                     dt,
                                                     localMatrix,
                                                     localRhs );
      }
    }
    else
    {
      // Take BCs defined for "pressure" field and apply values to "facePressure"
      fsManager.apply< FaceManager >( time_n + dt,
                                      mesh,
                                      fields::flow::pressure::key(),
                                      [&] ( FieldSpecificationBase const & fs,
                                            string const & setName,
                                            SortedArrayView< localIndex const > const & targetSet,
                                            FaceManager & targetGroup,
                                            string const & )
      {
        BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );

        if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
        {
          globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( stencil.size() );
          GEOS_LOG_RANK_0( GEOS_FMT( faceBcLogMessage,
                                     this->getName(), time_n+dt, fs.getCatalogName(), fs.getName(),
                                     setName, targetGroup.getName(), numTargetFaces ) );
        }

        if( stencil.size() == 0 )
        {
          return;
        }

        // first, evaluate BC to get primary field values (pressure)
        fs.applyFieldValue< FieldSpecificationEqual,
                            parallelDevicePolicy<> >( targetSet,
                                                      time_n + dt,
                                                      targetGroup,
                                                      fields::flow::facePressure::key() );


        // TODO: currently we just use model from the first cell in this stencil
        //       since it's not clear how to create fluid kernel wrappers for arbitrary models.
        //       Can we just use cell properties for an approximate flux computation?
        //       Then we can forget about capturing the fluid model.
        localIndex const er = stencil.getElementRegionIndices()( 0, 0 );
        localIndex const esr = stencil.getElementSubRegionIndices()( 0, 0 );
        ElementSubRegionBase & subRegion = mesh.getElemManager().getRegion( er ).getSubRegion( esr );
        string const & fluidName = subRegion.getReference< string >( BASE::viewKeyStruct::fluidNamesString() );
        SingleFluidBase & fluidBase = subRegion.getConstitutiveModel< SingleFluidBase >( fluidName );

        BoundaryStencilWrapper const stencilWrapper = stencil.createKernelWrapper();

        singlePhaseFVMKernels::
          DirichletFluxComputeKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                     dofKey,
                                                     this->getName(),
                                                     faceManager,
                                                     elemManager,
                                                     stencilWrapper,
                                                     fluidBase,
                                                     dt,
                                                     localMatrix,
                                                     localRhs );

      } );
    }
  } );
}

template<>
void SinglePhaseFVM< SinglePhaseProppantBase >::applyAquiferBC( real64 const GEOS_UNUSED_PARAM( time ),
                                                                real64 const GEOS_UNUSED_PARAM( dt ),
                                                                DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                                                DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                                CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                                                arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) ) const
{
  // Aquifer does not make sense for proppant flow in fractures
}

template<>
void SinglePhaseFVM<>::applyAquiferBC( real64 const time,
                                       real64 const dt,
                                       DomainPartition & domain,
                                       DofManager const & dofManager,
                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                       arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( this->getName() + "/accessors/" + elemDofKey );

    typename FluxComputeKernelBase::SinglePhaseFlowAccessors flowAccessors( elemManager, this->getName() );
    typename FluxComputeKernelBase::SinglePhaseFluidAccessors fluidAccessors( elemManager, this->getName() );

    fsManager.apply< FaceManager,
                     AquiferBoundaryCondition >( time + dt,
                                                 mesh,
                                                 AquiferBoundaryCondition::catalogName(),
                                                 [&] ( AquiferBoundaryCondition const & bc,
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

      AquiferBoundaryCondition::KernelWrapper aquiferBCWrapper = bc.createKernelWrapper();
      real64 const & aquiferDens = bc.getWaterPhaseDensity();

      singlePhaseFVMKernels::AquiferBCKernel::launch( stencil,
                                                      dofManager.rankOffset(),
                                                      elemDofNumber.toNestedViewConst(),
                                                      flowAccessors.get< fields::ghostRank >(),
                                                      aquiferBCWrapper,
                                                      aquiferDens,
                                                      flowAccessors.get< fields::flow::pressure >(),
                                                      flowAccessors.get< fields::flow::pressure_n >(),
                                                      flowAccessors.get< fields::flow::gravityCoefficient >(),
                                                      fluidAccessors.get< fields::singlefluid::density >(),
                                                      fluidAccessors.get< fields::singlefluid::dDensity_dPressure >(),
                                                      time,
                                                      dt,
                                                      localMatrix.toViewConstSizes(),
                                                      localRhs.toView() );

    } );
  } );

}

namespace
{
typedef SinglePhaseFVM< SinglePhaseProppantBase > SinglePhaseFVMProppant;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseFVMProppant, string const &, Group * const )
typedef SinglePhaseFVM<> SinglePhaseFVM;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseFVM, string const &, Group * const )
}
} /* namespace geos */
