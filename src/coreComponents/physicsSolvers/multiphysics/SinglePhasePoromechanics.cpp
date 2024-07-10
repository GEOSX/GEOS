/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhasePoromechanics.cpp
 */

#define GEOSX_DISPATCH_VEM /// enables VEM in FiniteElementDispatch

#include "SinglePhasePoromechanics.hpp"

#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/solid/PorousDamageSolid.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsDamage.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanics.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp"
#include "physicsSolvers/contact/SolidMechanicsLagrangeContact.hpp"
#include "physicsSolvers/contact/SolidMechanicsEmbeddedFractures.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::SinglePhasePoromechanics( const string & name,
                                                                                     Group * const parent )
  : Base( name, parent ),
  m_damageFlag()
{
  LinearSolverParameters & linearSolverParameters = this->m_linearSolverParameters.get();
  linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics;
  linearSolverParameters.mgr.separateComponents = true;
  linearSolverParameters.mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
  linearSolverParameters.dofsPerNode = 3;

  this->registerWrapper( viewKeyStruct::damageFlagString(), &m_damageFlag ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The flag to indicate whether a damage solid model is used" );
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::postInputInitialization()
{
  Base::postInputInitialization();

  GEOS_ERROR_IF( this->flowSolver()->getCatalogName() == "SinglePhaseReservoir" &&
                 this->getNonlinearSolverParameters().couplingType() != NonlinearSolverParameters::CouplingType::Sequential,
                 GEOS_FMT( "{}: {} solver is only designed to work for {} = {}",
                           this->getName(), catalogName(), NonlinearSolverParameters::viewKeysStruct::couplingTypeString(),
                           EnumStrings< NonlinearSolverParameters::CouplingType >::toString( NonlinearSolverParameters::CouplingType::Sequential )
                           ));
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::setupCoupling( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                                               DofManager & dofManager ) const
{
  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::setupSystem( DomainPartition & domain,
                                                                             DofManager & dofManager,
                                                                             CRSMatrix< real64, globalIndex > & localMatrix,
                                                                             ParallelVector & rhs,
                                                                             ParallelVector & solution,
                                                                             bool const setSparsity )
{
  if( this->m_precond )
  {
    this->m_precond->clear();
  }

  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  if( !this->m_precond && this->m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    createPreconditioner();
  }
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::initializePostInitialConditionsPreSubGroups()
{
  Base::initializePostInitialConditionsPreSubGroups();

  arrayView1d< string const > const & poromechanicsTargetRegionNames =
    this->template getReference< array1d< string > >( SolverBase::viewKeyStruct::targetRegionsString() );
  arrayView1d< string const > const & flowTargetRegionNames =
    this->flowSolver()->template getReference< array1d< string > >( SolverBase::viewKeyStruct::targetRegionsString() );
  for( integer i = 0; i < poromechanicsTargetRegionNames.size(); ++i )
  {
    GEOS_THROW_IF( std::find( flowTargetRegionNames.begin(), flowTargetRegionNames.end(), poromechanicsTargetRegionNames[i] )
                   == flowTargetRegionNames.end(),
                   GEOS_FMT( "{} {}: region `{}` must be a target region of `{}`",
                             getCatalogName(), this->getDataContext(), poromechanicsTargetRegionNames[i], this->flowSolver()->getDataContext() ),
                   InputError );
  }

  if( this->m_isThermal )
  {
    this->m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::thermalSinglePhasePoromechanics;
  }
  else
  {
    if( this->flowSolver()->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM )
    {
      this->m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::hybridSinglePhasePoromechanics;
    }
  }
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::assembleSystem( real64 const time_n,
                                                                                real64 const dt,
                                                                                DomainPartition & domain,
                                                                                DofManager const & dofManager,
                                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // Steps 1 and 2: compute element-based terms (mechanics and local flow terms)
  assembleElementBasedTerms( time_n,
                             dt,
                             domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  // Step 3: compute the fluxes (face-based contributions)
  if( m_stabilizationType == StabilizationType::Global || m_stabilizationType == StabilizationType::Local )
  {
    this->flowSolver()->assembleStabilizedFluxTerms( dt,
                                                     domain,
                                                     dofManager,
                                                     localMatrix,
                                                     localRhs );
  }
  else
  {
    this->flowSolver()->assembleFluxTerms( dt,
                                           domain,
                                           dofManager,
                                           localMatrix,
                                           localRhs );
  }

}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::assembleElementBasedTerms( real64 const time_n,
                                                                                           real64 const dt,
                                                                                           DomainPartition & domain,
                                                                                           DofManager const & dofManager,
                                                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                           arrayView1d< real64 > const & localRhs )
{
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( dt );

  real64 poromechanicsMaxForce = 0.0;
  real64 mechanicsMaxForce = 0.0;

  // step 1: apply the full poromechanics coupling on the target regions on the poromechanics solver

  set< string > poromechanicsRegionNames;

  this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                               MeshLevel & mesh,
                                                                               arrayView1d< string const > const & regionNames )
  {
    poromechanicsRegionNames.insert( regionNames.begin(), regionNames.end() );

    string const flowDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    if( m_damageFlag )
    {
      poromechanicsMaxForce =
        assemblyLaunch< constitutive::PorousDamageSolidBase,
                        poromechanicsDamageKernels::SinglePhasePoromechanicsDamageKernelFactory >( mesh,
                                                                                                   dofManager,
                                                                                                   regionNames,
                                                                                                   viewKeyStruct::porousMaterialNamesString(),
                                                                                                   localMatrix,
                                                                                                   localRhs,
                                                                                                   dt,
                                                                                                   flowDofKey,
                                                                                                   FlowSolverBase::viewKeyStruct::fluidNamesString() );
    }
    else if( this->m_isThermal )
    {
      poromechanicsMaxForce =
        assemblyLaunch< constitutive::PorousSolid< ElasticIsotropic >, // TODO: change once there is a cmake solution
                        thermalPoromechanicsKernels::ThermalSinglePhasePoromechanicsKernelFactory >( mesh,
                                                                                                     dofManager,
                                                                                                     regionNames,
                                                                                                     viewKeyStruct::porousMaterialNamesString(),
                                                                                                     localMatrix,
                                                                                                     localRhs,
                                                                                                     dt,
                                                                                                     flowDofKey,
                                                                                                     this->m_performStressInitialization,
                                                                                                     FlowSolverBase::viewKeyStruct::fluidNamesString() );
    }
    else
    {
      poromechanicsMaxForce =
        assemblyLaunch< constitutive::PorousSolidBase,
                        poromechanicsKernels::SinglePhasePoromechanicsKernelFactory >( mesh,
                                                                                       dofManager,
                                                                                       regionNames,
                                                                                       viewKeyStruct::porousMaterialNamesString(),
                                                                                       localMatrix,
                                                                                       localRhs,
                                                                                       dt,
                                                                                       flowDofKey,
                                                                                       this->m_performStressInitialization,
                                                                                       FlowSolverBase::viewKeyStruct::fluidNamesString() );
    }
  } );

  // step 2: apply mechanics solver on its target regions not included in the poromechanics solver target regions

  this->solidMechanicsSolver()->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                                              MeshLevel & mesh,
                                                                                              arrayView1d< string const > const & regionNames )
  {
    // collect the target region of the mechanics solver not included in the poromechanics target regions
    array1d< string > filteredRegionNames;
    filteredRegionNames.reserve( regionNames.size() );
    for( string const & regionName : regionNames )
    {
      // if the mechanics target region is not included in the poromechanics target region, save the string
      if( poromechanicsRegionNames.count( regionName ) == 0 )
      {
        filteredRegionNames.emplace_back( regionName );
      }
    }

    // if the array is empty, the mechanics and poromechanics solver target regions overlap perfectly, there is nothing to do
    if( filteredRegionNames.empty() )
    {
      return;
    }

    mechanicsMaxForce =
      assemblyLaunch< constitutive::SolidBase,
                      solidMechanicsLagrangianFEMKernels::QuasiStaticFactory >( mesh,
                                                                                dofManager,
                                                                                filteredRegionNames.toViewConst(),
                                                                                SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString(),
                                                                                localMatrix,
                                                                                localRhs,
                                                                                dt );
  } );

  this->solidMechanicsSolver()->applyContactConstraint( dofManager, domain, localMatrix, localRhs );
  this->solidMechanicsSolver()->getMaxForce() = LvArray::math::max( mechanicsMaxForce, poromechanicsMaxForce );
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::createPreconditioner()
{
  if( this->m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::UpperTriangular,
                                                                           SchurComplementOption::RowsumDiagonalProbing,
                                                                           BlockScalingOption::FrobeniusNorm );

    auto mechPrecond = LAInterface::createPreconditioner( this->solidMechanicsSolver()->getLinearSolverParameters() );
    precond->setupBlock( 0,
                         { { solidMechanics::totalDisplacement::key(), { 3, true } } },
                         std::make_unique< SeparateComponentPreconditioner< LAInterface > >( 3, std::move( mechPrecond ) ) );

    auto flowPrecond = LAInterface::createPreconditioner( this->flowSolver()->getLinearSolverParameters() );
    precond->setupBlock( 1,
                         { { flow::pressure::key(), { 1, true } } },
                         std::move( flowPrecond ) );

    this->m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                               MeshLevel & mesh,
                                                                               arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   CellElementSubRegion & subRegion )
    {
      this->flowSolver()->updateFluidState( subRegion );
      if( this->m_isThermal )
      {
        this->flowSolver()->updateSolidInternalEnergyModel( subRegion );
      }
    } );
  } );
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void SinglePhasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::updateBulkDensity( ElementSubRegionBase & subRegion )
{
  // get the fluid model (to access fluid density)
  string const fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
  SingleFluidBase const & fluid = this->template getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

  // get the solid model (to access porosity and solid density)
  string const solidName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
  CoupledSolidBase const & solid = this->template getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

  // update the bulk density
  poromechanicsKernels::
    SinglePhaseBulkDensityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( fluid,
                                               solid,
                                               subRegion );
}

template class SinglePhasePoromechanics<>;
template class SinglePhasePoromechanics< SinglePhaseBase, SolidMechanicsLagrangeContact >;
template class SinglePhasePoromechanics< SinglePhaseBase, SolidMechanicsEmbeddedFractures >;
template class SinglePhasePoromechanics< SinglePhaseReservoirAndWells<> >;
//template class SinglePhasePoromechanics< SinglePhaseReservoirAndWells<>, SolidMechanicsLagrangeContact >;
//template class SinglePhasePoromechanics< SinglePhaseReservoirAndWells<>, SolidMechanicsEmbeddedFractures >;

namespace
{
typedef SinglePhasePoromechanics< SinglePhaseReservoirAndWells<> > SinglePhaseReservoirPoromechanics;
REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseReservoirPoromechanics, string const &, Group * const )
typedef SinglePhasePoromechanics<> SinglePhasePoromechanics;
REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanics, string const &, Group * const )
}

} /* namespace geos */
