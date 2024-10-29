/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
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
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "mesh/PerforationFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
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
{}

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
  if( flowSolver()->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseHybridFVM )
  {
    // add Reservoir
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseReservoirHybridFVM;
  }
  else
  {
    // add Reservoir
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseReservoirFVM;
  }
}

template<>
void
CompositionalMultiphaseReservoirAndWells< MultiphasePoromechanics<> >::
setMGRStrategy()
{
  // flow solver here is indeed flow solver, not poromechanics solver
  if( flowSolver()->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseHybridFVM )
  {
    GEOS_LOG_RANK_0( "The poromechanics MGR strategy for hybrid FVM is not implemented" );
  }
  else
  {
    // add Reservoir
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::multiphasePoromechanicsReservoirFVM;
  }
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
initializePostInitialConditionsPreSubGroups()
{
  Base::initializePostInitialConditionsPreSubGroups();
  setMGRStrategy();
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

    // TODO: remove this and just call SolverBase::setupSystem when DofManager can handle the coupling

    // Populate off-diagonal sparsity between well and reservoir

    integer const resNDOF = Base::wellSolver()->numDofPerResElement();
    integer const wellNDOF = Base::wellSolver()->numDofPerWellElement();

    integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
    integer constexpr maxNumDof  = maxNumComp + 1;

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

  using TAG = compositionalMultiphaseWellKernels::SubRegionTag;
  using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = compositionalMultiphaseWellKernels::ColOffset;

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

    integer constexpr MAX_NUM_COMP = MultiFluidBase::MAX_NUM_COMPONENTS;
    integer constexpr MAX_NUM_DOF = MAX_NUM_COMP + 1;

    integer const numComps = Base::wellSolver()->numFluidComponents();
    integer const resNumDofs = Base::wellSolver()->numDofPerResElement();

    string const resDofKey = dofManager.getKey( Base::wellSolver()->resElementDofName() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const resDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );
    ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber =
      resDofNumberAccessor.toNestedViewConst();
    globalIndex const rankOffset = dofManager.rankOffset();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion const & subRegion )
    {

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

      areWellsShut = 0;

      PerforationData const * const perforationData = subRegion.getPerforationData();

      // get the degrees of freedom
      string const wellDofKey = dofManager.getKey( Base::wellSolver()->wellElementDofName() );
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      // get well variables on perforations
      arrayView2d< real64 const > const & compPerfRate =
        perforationData->getField< fields::well::compPerforationRate >();
      arrayView3d< real64 const > const & dCompPerfRate_dPres =
        perforationData->getField< fields::well::dCompPerforationRate_dPres >();
      arrayView4d< real64 const > const & dCompPerfRate_dComp =
        perforationData->getField< fields::well::dCompPerforationRate_dComp >();

      arrayView1d< localIndex const > const & perfWellElemIndex =
        perforationData->getField< fields::perforation::wellElementIndex >();

      // get the element region, subregion, index
      arrayView1d< localIndex const > const & resElementRegion =
        perforationData->getField< fields::perforation::reservoirElementRegion >();
      arrayView1d< localIndex const > const & resElementSubRegion =
        perforationData->getField< fields::perforation::reservoirElementSubRegion >();
      arrayView1d< localIndex const > const & resElementIndex =
        perforationData->getField< fields::perforation::reservoirElementIndex >();

      bool const useTotalMassEquation = this->flowSolver()->useTotalMassEquation() > 0;

      RAJA::ReduceSum< parallelDeviceReduce, integer > numCrossflowPerforations( 0 );

      // loop over the perforations and add the rates to the residual and jacobian
      forAll< parallelDevicePolicy<> >( perforationData->size(), [=] GEOS_HOST_DEVICE ( localIndex const iperf )
      {
        // local working variables and arrays
        stackArray1d< localIndex, 2 * MAX_NUM_COMP > eqnRowIndices( 2 * numComps );
        stackArray1d< globalIndex, 2 * MAX_NUM_DOF > dofColIndices( 2 * resNumDofs );

        stackArray1d< real64, 2 * MAX_NUM_COMP > localPerf( 2 * numComps );
        stackArray2d< real64, 2 * MAX_NUM_COMP * 2 * MAX_NUM_DOF > localPerfJacobian( 2 * numComps, 2 * resNumDofs );

        // get the reservoir (sub)region and element indices
        localIndex const er  = resElementRegion[iperf];
        localIndex const esr = resElementSubRegion[iperf];
        localIndex const ei  = resElementIndex[iperf];

        // get the well element index for this perforation
        localIndex const iwelem = perfWellElemIndex[iperf];
        globalIndex const resOffset = resDofNumber[er][esr][ei];
        globalIndex const wellElemOffset = wellElemDofNumber[iwelem];

        for( integer ic = 0; ic < numComps; ++ic )
        {
          eqnRowIndices[TAG::RES * numComps + ic] = LvArray::integerConversion< localIndex >( resOffset - rankOffset ) + ic;
          eqnRowIndices[TAG::WELL * numComps + ic] = LvArray::integerConversion< localIndex >( wellElemOffset - rankOffset ) + ROFFSET::MASSBAL + ic;
        }
        for( integer jdof = 0; jdof < resNumDofs; ++jdof )
        {
          dofColIndices[TAG::RES * resNumDofs + jdof] = resOffset + jdof;
          dofColIndices[TAG::WELL * resNumDofs + jdof] = wellElemOffset + COFFSET::DPRES + jdof;
        }

        // populate local flux vector and derivatives
        for( integer ic = 0; ic < numComps; ++ic )
        {
          localPerf[TAG::RES * numComps + ic] = dt * compPerfRate[iperf][ic];
          localPerf[TAG::WELL * numComps + ic] = -dt * compPerfRate[iperf][ic];

          if( detectCrossflow )
          {
            if( compPerfRate[iperf][ic] > LvArray::NumericLimits< real64 >::epsilon )
            {
              numCrossflowPerforations += 1;
            }
          }

          for( integer ke = 0; ke < 2; ++ke )
          {
            localIndex const localDofIndexPres = ke * resNumDofs;
            localPerfJacobian[TAG::RES * numComps + ic][localDofIndexPres] = dt * dCompPerfRate_dPres[iperf][ke][ic];
            localPerfJacobian[TAG::WELL * numComps + ic][localDofIndexPres] = -dt * dCompPerfRate_dPres[iperf][ke][ic];

            for( integer jc = 0; jc < numComps; ++jc )
            {
              localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
              localPerfJacobian[TAG::RES * numComps + ic][localDofIndexComp] = dt * dCompPerfRate_dComp[iperf][ke][ic][jc];
              localPerfJacobian[TAG::WELL * numComps + ic][localDofIndexComp] = -dt * dCompPerfRate_dComp[iperf][ke][ic][jc];
            }
          }
        }

        if( useTotalMassEquation )
        {
          // Apply equation/variable change transformation(s)
          stackArray1d< real64, 2 * MAX_NUM_DOF > work( 2 * resNumDofs );
          shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComps, numComps, resNumDofs * 2, 2, localPerfJacobian, work );
          shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( numComps, numComps, 2, localPerf );
        }

        for( localIndex i = 0; i < localPerf.size(); ++i )
        {
          if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < localMatrix.numRows() )
          {
            localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                              dofColIndices.data(),
                                                                              localPerfJacobian[i].dataIfContiguous(),
                                                                              2 * resNumDofs );
            RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndices[i]], localPerf[i] );
          }
        }
      } );


      if( detectCrossflow ) // check to avoid communications if not needed
      {
        globalIndex const totalNumCrossflowPerforations = MpiWrapper::sum( numCrossflowPerforations.get() );
        if( totalNumCrossflowPerforations > 0 )
        {
          GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "CompositionalMultiphaseReservoir '{}': Warning! Crossflow detected at {} perforations in well {}"
                                              "To disable crossflow for injectors, you can use the field '{}' in the WellControls '{}' section",
                                              this->getName(), totalNumCrossflowPerforations, subRegion.getName(),
                                              WellControls::viewKeyStruct::enableCrossflowString(), wellControls.getName() ) );
        }
      }
    } );

    // update dynamically the MGR recipe to optimize the linear solve if all wells are shut
    areWellsShut = MpiWrapper::min( areWellsShut );
    m_linearSolverParameters.get().mgr.areWellsShut = areWellsShut;

  } );
}

template class CompositionalMultiphaseReservoirAndWells<>;
template class CompositionalMultiphaseReservoirAndWells< MultiphasePoromechanics<> >;

namespace
{
typedef CompositionalMultiphaseReservoirAndWells<> CompositionalMultiphaseFlowAndWells;
typedef CompositionalMultiphaseReservoirAndWells< MultiphasePoromechanics<> > CompositionalMultiphasePoromechanicsAndWells;
REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseFlowAndWells, string const &, Group * const )
REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphasePoromechanicsAndWells, string const &, Group * const )
}

} /* namespace geos */
