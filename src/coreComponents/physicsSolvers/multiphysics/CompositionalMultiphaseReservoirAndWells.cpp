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
 * @file CompositionalMultiphaseReservoirAndWells.cpp
 *
 */

#include "CompositionalMultiphaseReservoirAndWells.hpp"

#include "common/TimingMacros.hpp"
#include "dataRepository/LogLevelsInfo.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "mesh/PerforationFields.hpp"
#include "physicsSolvers/multiphysics/CoupledReservoirAndWellKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#include "physicsSolvers/multiphysics/MultiphasePoromechanics.hpp"


namespace geos
{

using namespace dataRepository;
using namespace constitutive;

template< typename RESERVOIR_SOLVER >
CompositionalMultiphaseReservoirAndWells< RESERVOIR_SOLVER >::
CompositionalMultiphaseReservoirAndWells( const string & name,
                                          Group * const parent )
  : Base( name, parent )
{
  Base::template addLogLevel< logInfo::Crossflow >();
}

template< typename RESERVOIR_SOLVER >
CompositionalMultiphaseReservoirAndWells< RESERVOIR_SOLVER >::
~CompositionalMultiphaseReservoirAndWells()
{}

template<>
CompositionalMultiphaseBase *
CompositionalMultiphaseReservoirAndWells<>::
flowSolver() const
{
  return this->reservoirSolver();
}

template<>
CompositionalMultiphaseBase *
CompositionalMultiphaseReservoirAndWells< MultiphasePoromechanics<> >::
flowSolver() const
{
  return this->reservoirSolver()->flowSolver();
}

template<>
void
CompositionalMultiphaseReservoirAndWells<>::
setMGRStrategy()
{
  LinearSolverParameters & linearSolverParameters = m_linearSolverParameters.get();

  if( linearSolverParameters.preconditionerType != LinearSolverParameters::PreconditionerType::mgr )
    return;

  linearSolverParameters.mgr.separateComponents = true;
  linearSolverParameters.dofsPerNode = 3;

  if( dynamic_cast< CompositionalMultiphaseHybridFVM * >( this->flowSolver() ) )
  {
    if( isThermal() )
    {
      GEOS_ERROR( GEOS_FMT( "{}: MGR strategy is not implemented for thermal {}/{}",
      this->getName(), this->getCatalogName(), this->flowSolver()->getCatalogName()));
    }
    else
    {
      linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseReservoirHybridFVM;
    }
  }
  else
  {
    if( isThermal() )
    {
      m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::thermalCompositionalMultiphaseReservoirFVM;
    }
    else
    {
      linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseReservoirFVM;
    }
  }
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}: MGR strategy set to {}", getName(),
                                      EnumStrings< LinearSolverParameters::MGR::StrategyType >::toString( linearSolverParameters.mgr.strategy )));
}

template<>
void
CompositionalMultiphaseReservoirAndWells< MultiphasePoromechanics<> >::
setMGRStrategy()
{
  LinearSolverParameters & linearSolverParameters = m_linearSolverParameters.get();

  if( linearSolverParameters.preconditionerType != LinearSolverParameters::PreconditionerType::mgr )
    return;

  linearSolverParameters.mgr.separateComponents = true;
  linearSolverParameters.dofsPerNode = 3;

  if( dynamic_cast< CompositionalMultiphaseHybridFVM * >( this->flowSolver() ) )
  {
    GEOS_ERROR( GEOS_FMT( "{}: MGR strategy is not implemented for {}/{}",
    this->getName(), this->getCatalogName(), this->flowSolver()->getCatalogName() ) );
  }
  else
  {
    linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::multiphasePoromechanicsReservoirFVM;
  }
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}: MGR strategy set to {}", getName(),
                                      EnumStrings< LinearSolverParameters::MGR::StrategyType >::toString( linearSolverParameters.mgr.strategy )));
}

template< typename RESERVOIR_SOLVER >
void
CompositionalMultiphaseReservoirAndWells< RESERVOIR_SOLVER >::
initializePreSubGroups()
{
  Base::initializePreSubGroups();

  CompositionalMultiphaseBase const * const flowSolver = this->flowSolver();
  Base::wellSolver()->setFlowSolverName( flowSolver->getName() );

  bool const useMassFlow = flowSolver->getReference< integer >( CompositionalMultiphaseBase::viewKeyStruct::useMassFlagString() );;
  bool const useMassWell = Base::wellSolver()->template getReference< integer >( CompositionalMultiphaseWell::viewKeyStruct::useMassFlagString() );
  GEOS_THROW_IF( useMassFlow != useMassWell,
                 GEOS_FMT( "{}: the input flag {} must be the same in the flow and well solvers, respectively '{}' and '{}'",
                           this->getDataContext(), CompositionalMultiphaseBase::viewKeyStruct::useMassFlagString(),
                           Base::reservoirSolver()->getDataContext(), Base::wellSolver()->getDataContext() ),
                 InputError );
}

template< typename RESERVOIR_SOLVER >
void
CompositionalMultiphaseReservoirAndWells< RESERVOIR_SOLVER >::
addCouplingSparsityPattern( DomainPartition const & domain,
                            DofManager const & dofManager,
                            SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOS_MARK_FUNCTION;

  this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                               MeshLevel const & mesh,
                                                                               arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    // TODO: remove this and just call PhysicsSolverBase::setupSystem when DofManager can handle the coupling

    // Populate off-diagonal sparsity between well and reservoir

    integer const resNDOF = Base::wellSolver()->numDofPerResElement();
    integer const wellNDOF = Base::wellSolver()->numDofPerWellElement();

    integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
    integer constexpr maxNumDof  = maxNumComp + 2;

    string const wellDofKey = dofManager.getKey( Base::wellSolver()->wellElementDofName() );
    string const resDofKey  = dofManager.getKey( Base::wellSolver()->resElementDofName() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const & resDofNumber =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );

    globalIndex const rankOffset = dofManager.rankOffset();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion const & subRegion )
    {
      PerforationData const * const perforationData = subRegion.getPerforationData();

      // get the well degrees of freedom and ghosting info
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      // get the well element indices corresponding to each perforation
      arrayView1d< localIndex const > const & perfWellElemIndex =
        perforationData->getField< fields::perforation::wellElementIndex >();

      // get the element region, subregion, index
      arrayView1d< localIndex const > const & resElementRegion =
        perforationData->getField< fields::perforation::reservoirElementRegion >();
      arrayView1d< localIndex const > const & resElementSubRegion =
        perforationData->getField< fields::perforation::reservoirElementSubRegion >();
      arrayView1d< localIndex const > const & resElementIndex =
        perforationData->getField< fields::perforation::reservoirElementIndex >();

      // Insert the entries corresponding to reservoir-well perforations
      // This will fill J_WR, and J_RW
      forAll< serialPolicy >( perforationData->size(), [=] ( localIndex const iperf )
      {
        stackArray1d< globalIndex, maxNumDof > eqnRowIndicesRes( resNDOF );
        stackArray1d< globalIndex, maxNumDof > eqnRowIndicesWell( wellNDOF );
        stackArray1d< globalIndex, maxNumDof > dofColIndicesRes( resNDOF );
        stackArray1d< globalIndex, maxNumDof > dofColIndicesWell( wellNDOF );

        // get the reservoir (sub)region and element indices
        localIndex const er = resElementRegion[iperf];
        localIndex const esr = resElementSubRegion[iperf];
        localIndex const ei = resElementIndex[iperf];
        localIndex const iwelem = perfWellElemIndex[iperf];

        for( integer idof = 0; idof < resNDOF; ++idof )
        {
          eqnRowIndicesRes[idof] = resDofNumber[er][esr][ei] + idof - rankOffset;
          dofColIndicesRes[idof] = resDofNumber[er][esr][ei] + idof;
        }

        for( integer idof = 0; idof < wellNDOF; ++idof )
        {
          eqnRowIndicesWell[idof] = wellElemDofNumber[iwelem] + idof - rankOffset;
          dofColIndicesWell[idof] = wellElemDofNumber[iwelem] + idof;
        }

        for( localIndex i = 0; i < eqnRowIndicesRes.size(); ++i )
        {
          if( eqnRowIndicesRes[i] >= 0 && eqnRowIndicesRes[i] < pattern.numRows() )
          {
            for( localIndex j = 0; j < dofColIndicesWell.size(); ++j )
            {
              pattern.insertNonZero( eqnRowIndicesRes[i], dofColIndicesWell[j] );
            }
          }
        }

        for( localIndex i = 0; i < eqnRowIndicesWell.size(); ++i )
        {
          if( eqnRowIndicesWell[i] >= 0 && eqnRowIndicesWell[i] < pattern.numRows() )
          {
            for( localIndex j = 0; j < dofColIndicesRes.size(); ++j )
            {
              pattern.insertNonZero( eqnRowIndicesWell[i], dofColIndicesRes[j] );
            }
          }
        }
      } );
    } );
  } );
}

template< typename RESERVOIR_SOLVER >
void
CompositionalMultiphaseReservoirAndWells< RESERVOIR_SOLVER >::
assembleCouplingTerms( real64 const time_n,
                       real64 const dt,
                       DomainPartition const & domain,
                       DofManager const & dofManager,
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs )
{
  using namespace compositionalMultiphaseUtilities;

  GEOS_THROW_IF( !Base::m_isWellTransmissibilityComputed,
                 GEOS_FMT( "{} {}: The well transmissibility has not been computed yet",
                           this->getCatalogName(), this->getName() ),
                 std::runtime_error );

  this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                               MeshLevel const & mesh,
                                                                               arrayView1d< string const > const & regionNames )
  {
    integer areWellsShut = 1;

    ElementRegionManager const & elemManager = mesh.getElemManager();

    integer const numComps = Base::wellSolver()->numFluidComponents();

    string const resDofKey = dofManager.getKey( Base::wellSolver()->resElementDofName() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const resDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );
    ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber =
      resDofNumberAccessor.toNestedViewConst();
    globalIndex const rankOffset = dofManager.rankOffset();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion const & subRegion )
    {
      string const & fluidName = this->flowSolver()->template getConstitutiveName< MultiFluidBase >( subRegion );
      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

      // if the well is shut, we neglect reservoir-well flow that may occur despite the zero rate
      // therefore, we do not want to compute perforation rates and we simply assume they are zero
      WellControls const & wellControls = Base::wellSolver()->getWellControls( subRegion );
      bool const detectCrossflow =
        ( wellControls.isInjector() ) && wellControls.isCrossflowEnabled() &&
        getLogLevel() >= 1; // since detect crossflow requires communication, we detect it only if the logLevel is sufficiently high

      if( !wellControls.isWellOpen( time_n + dt ) )
      {
        return;
      }

      PerforationData const * const perforationData = subRegion.getPerforationData();

      // get the degrees of freedom
      string const wellDofKey = dofManager.getKey( Base::wellSolver()->wellElementDofName() );
      areWellsShut = 0;

      integer useTotalMassEquation=Base::wellSolver()->useTotalMassEquation();
      integer numCrossflowPerforations=0;
      if( isThermal ( )  )
      {
        coupledReservoirAndWellKernels::
          ThermalCompositionalMultiPhaseFluxKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( numComps,
                                                     wellControls.isProducer(),
                                                     dt,
                                                     rankOffset,
                                                     wellDofKey,
                                                     subRegion,
                                                     resDofNumber,
                                                     perforationData,
                                                     fluid,
                                                     useTotalMassEquation,
                                                     detectCrossflow,
                                                     numCrossflowPerforations,
                                                     localRhs,
                                                     localMatrix );
      }
      else
      {
        coupledReservoirAndWellKernels::
          IsothermalCompositionalMultiPhaseFluxKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( numComps,
                                                     dt,
                                                     rankOffset,
                                                     wellDofKey,
                                                     subRegion,
                                                     resDofNumber,
                                                     perforationData,
                                                     fluid,
                                                     useTotalMassEquation,
                                                     detectCrossflow,
                                                     numCrossflowPerforations,
                                                     localRhs,
                                                     localMatrix );
      }

      if( detectCrossflow )                                                         // check to avoid communications if not needed
      {
        globalIndex const totalNumCrossflowPerforations = MpiWrapper::sum( numCrossflowPerforations );
        if( totalNumCrossflowPerforations > 0 )
        {
          GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Crossflow, GEOS_FMT( "CompositionalMultiphaseReservoir '{}': Warning! Crossflow detected at {} perforations in well {}"
                                                                    "To disable crossflow for injectors, you can use the field '{}' in the WellControls '{}' section",
                                                                    this->getName(), totalNumCrossflowPerforations, subRegion.getName(),
                                                                    WellControls::viewKeyStruct::enableCrossflowString(), wellControls.getName() ));
        }
      }


      // update dynamically the MGR recipe to optimize the linear solve if all wells are shut
      areWellsShut = MpiWrapper::min( areWellsShut );
      m_linearSolverParameters.get().mgr.areWellsShut = areWellsShut;
    } );
  } );
}

template class CompositionalMultiphaseReservoirAndWells<>;
template class CompositionalMultiphaseReservoirAndWells< MultiphasePoromechanics<> >;

namespace
{
typedef CompositionalMultiphaseReservoirAndWells<> CompositionalMultiphaseFlowAndWells;
typedef CompositionalMultiphaseReservoirAndWells< MultiphasePoromechanics<> > CompositionalMultiphasePoromechanicsAndWells;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, CompositionalMultiphaseFlowAndWells, string const &, Group * const )
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, CompositionalMultiphasePoromechanicsAndWells, string const &, Group * const )
}

} /* namespace geos */
