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
 * @file MultiphasePoromechanics.cpp
 */

#define GEOS_DISPATCH_VEM /// enables VEM in FiniteElementDispatch

#include "MultiphasePoromechanics.hpp"

#include "dataRepository/LogLevelsInfo.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/MultiphasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalMultiphasePoromechanics.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp"
//#include "physicsSolvers/contact/SolidMechanicsLagrangeContact.hpp"
//#include "physicsSolvers/contact/SolidMechanicsEmbeddedFractures.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;
using namespace stabilization;

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
MultiphasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::MultiphasePoromechanics( const string & name,
                                                                                   Group * const parent )
  : Base( name, parent )
{}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::postInputInitialization()
{
  Base::postInputInitialization();

  setMGRStrategy();
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::setupCoupling( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                                              DofManager & dofManager ) const
{
  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::assembleSystem( real64 const time,
                                                                               real64 const dt,
                                                                               DomainPartition & domain,
                                                                               DofManager const & dofManager,
                                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                               arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // Steps 1 and 2: compute element-based terms (mechanics and local flow terms)
  assembleElementBasedTerms( time,
                             dt,
                             domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  // step 3: compute the fluxes (face-based contributions)

  if( m_stabilizationType == StabilizationType::Global ||
      m_stabilizationType == StabilizationType::Local )
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

template<>
void MultiphasePoromechanics< CompositionalMultiphaseReservoirAndWells<>, SolidMechanicsLagrangianFEM >::assembleSystem( real64 const time,
                                                                                                                         real64 const dt,
                                                                                                                         DomainPartition & domain,
                                                                                                                         DofManager const & dofManager,
                                                                                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                                                         arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // Steps 1 and 2: compute element-based terms (mechanics and local flow terms)
  assembleElementBasedTerms( time,
                             dt,
                             domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  // step 3: compute the fluxes (face-based contributions)

  if( m_stabilizationType == StabilizationType::Global ||
      m_stabilizationType == StabilizationType::Local )
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

  // step 4: assemble well contributions

  this->flowSolver()->wellSolver()->assembleSystem( time, dt, domain, dofManager, localMatrix, localRhs );
  this->flowSolver()->assembleCouplingTerms( time, dt, domain, dofManager, localMatrix, localRhs );
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::assembleElementBasedTerms( real64 const time_n,
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

    string const flowDofKey = dofManager.getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() );

    if( this->m_isThermal )
    {
      poromechanicsMaxForce =
        assemblyLaunch< constitutive::PorousSolidBase,
                        thermalPoromechanicsKernels::ThermalMultiphasePoromechanicsKernelFactory >( mesh,
                                                                                                    dofManager,
                                                                                                    regionNames,
                                                                                                    Base::viewKeyStruct::porousMaterialNamesString(),
                                                                                                    localMatrix,
                                                                                                    localRhs,
                                                                                                    dt,
                                                                                                    flowDofKey,
                                                                                                    this->flowSolver()->numFluidComponents(),
                                                                                                    this->flowSolver()->numFluidPhases(),
                                                                                                    this->flowSolver()->useTotalMassEquation(),
                                                                                                    this->m_performStressInitialization,
                                                                                                    FlowSolverBase::viewKeyStruct::fluidNamesString() );
    }
    else
    {
      poromechanicsMaxForce =
        assemblyLaunch< constitutive::PorousSolidBase,
                        poromechanicsKernels::MultiphasePoromechanicsKernelFactory >( mesh,
                                                                                      dofManager,
                                                                                      regionNames,
                                                                                      Base::viewKeyStruct::porousMaterialNamesString(),
                                                                                      localMatrix,
                                                                                      localRhs,
                                                                                      dt,
                                                                                      flowDofKey,
                                                                                      this->flowSolver()->numFluidComponents(),
                                                                                      this->flowSolver()->numFluidPhases(),
                                                                                      this->flowSolver()->useSimpleAccumulation(),
                                                                                      this->flowSolver()->useTotalMassEquation(),
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


  this->solidMechanicsSolver()->getMaxForce() = LvArray::math::max( mechanicsMaxForce, poromechanicsMaxForce );
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::initializePostInitialConditionsPreSubGroups()
{
  Base::initializePostInitialConditionsPreSubGroups();

  arrayView1d< string const > const & poromechanicsTargetRegionNames =
    this->template getReference< array1d< string > >( PhysicsSolverBase::viewKeyStruct::targetRegionsString() );
  arrayView1d< string const > const & solidMechanicsTargetRegionNames =
    this->solidMechanicsSolver()->template getReference< array1d< string > >( PhysicsSolverBase::viewKeyStruct::targetRegionsString() );
  arrayView1d< string const > const & flowTargetRegionNames =
    this->flowSolver()->template getReference< array1d< string > >( PhysicsSolverBase::viewKeyStruct::targetRegionsString() );
  for( integer i = 0; i < poromechanicsTargetRegionNames.size(); ++i )
  {
    GEOS_THROW_IF( std::find( solidMechanicsTargetRegionNames.begin(), solidMechanicsTargetRegionNames.end(),
                              poromechanicsTargetRegionNames[i] )
                   == solidMechanicsTargetRegionNames.end(),
                   GEOS_FMT( "{} {}: region {} must be a target region of {}",
                             getCatalogName(), this->getDataContext(), poromechanicsTargetRegionNames[i],
                             this->solidMechanicsSolver()->getDataContext() ),
                   InputError );
    GEOS_THROW_IF( std::find( flowTargetRegionNames.begin(), flowTargetRegionNames.end(), poromechanicsTargetRegionNames[i] )
                   == flowTargetRegionNames.end(),
                   GEOS_FMT( "{} {}: region `{}` must be a target region of `{}`",
                             getCatalogName(), this->getDataContext(), poromechanicsTargetRegionNames[i], this->flowSolver()->getDataContext() ),
                   InputError );
  }
}

template<>
void MultiphasePoromechanics<>::setMGRStrategy()
{
  LinearSolverParameters & linearSolverParameters = this->m_linearSolverParameters.get();

  if( linearSolverParameters.preconditionerType != LinearSolverParameters::PreconditionerType::mgr )
    return;

  linearSolverParameters.mgr.separateComponents = true;
  linearSolverParameters.dofsPerNode = 3;

  if( dynamic_cast< CompositionalMultiphaseHybridFVM * >( this->flowSolver() ) )
  {
    GEOS_ERROR( GEOS_FMT( "{}: MGR strategy is not implemented for {}/{}",
                          this->getName(), this->getCatalogName(), this->flowSolver()->getCatalogName()));
  }
  else
  {
    if( this->m_isThermal )
    {
      linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::thermalMultiphasePoromechanics;
    }
    else
    {
      linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::multiphasePoromechanics;
    }
  }
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}: MGR strategy set to {}", getName(),
                                      EnumStrings< LinearSolverParameters::MGR::StrategyType >::toString( linearSolverParameters.mgr.strategy )));
}

template<>
void MultiphasePoromechanics< CompositionalMultiphaseReservoirAndWells<>, SolidMechanicsLagrangianFEM >::setMGRStrategy()
{
  // same as CompositionalMultiphaseReservoirAndWells< MultiphasePoromechanics<> >::setMGRStrategy

  LinearSolverParameters & linearSolverParameters = m_linearSolverParameters.get();

  if( linearSolverParameters.preconditionerType != LinearSolverParameters::PreconditionerType::mgr )
    return;

  // flow solver here is indeed flow solver, not poromechanics solver
  if( dynamic_cast< CompositionalMultiphaseHybridFVM * >( flowSolver() ) )
  {
    GEOS_ERROR( "The poromechanics MGR strategy for hybrid FVM is not implemented" );
  }
  else
  {
    // add Reservoir
    linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::multiphasePoromechanicsReservoirFVM;
  }
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}: MGR strategy set to {}", getName(),
                                      EnumStrings< LinearSolverParameters::MGR::StrategyType >::toString( linearSolverParameters.mgr.strategy )));
  // TODO check if needed
  //linearSolverParameters.mgr.separateComponents = true;
  //linearSolverParameters.mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
}

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER, MECHANICS_SOLVER >::updateBulkDensity( ElementSubRegionBase & subRegion )
{
  // get the fluid model (to access fluid density)
  string const fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
  MultiFluidBase const & fluid = this->template getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

  // get the solid model (to access porosity and solid density)
  string const solidName = subRegion.getReference< string >( Base::viewKeyStruct::porousMaterialNamesString() );
  CoupledSolidBase const & solid = this->template getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

  // update the bulk density
  poromechanicsKernels::
    MultiphaseBulkDensityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( this->flowSolver()->numFluidComponents(),
                                               fluid,
                                               solid,
                                               subRegion );
}

template class MultiphasePoromechanics<>;
//template class MultiphasePoromechanics< CompositionalMultiphaseBase, SolidMechanicsLagrangeContact >;
//template class MultiphasePoromechanics< CompositionalMultiphaseBase, SolidMechanicsEmbeddedFractures >;
template class MultiphasePoromechanics< CompositionalMultiphaseReservoirAndWells<> >;

namespace
{
typedef MultiphasePoromechanics< CompositionalMultiphaseReservoirAndWells<> > MultiphaseReservoirPoromechanics;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, MultiphaseReservoirPoromechanics, string const &, Group * const )
typedef MultiphasePoromechanics<> MultiphasePoromechanics;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, MultiphasePoromechanics, string const &, Group * const )
}

} /* namespace geos */
