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
 * @file CoupledReservoirAndWellsBase.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDRESERVOIRANDWELLSBASE_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDRESERVOIRANDWELLSBASE_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBase.hpp"

namespace geos
{

namespace coupledReservoirAndWellsInternal
{
/**
 * @brief Utility function for the implementation details of addCouplingNumZeros
 * @param solver the coupled solver
 * @param domain the physical domain object
 * @param dofManager degree-of-freedom manager associated with the linear system
 * @param rowLengths the row-by-row length
 * @param resNumDof number of reservoir element dofs
 * @param wellNumDof number of well element dofs
 * @param resElemDofName name of the reservoir element dofs
 * @param wellElemDofName name of the well element dofs
 */
void
addCouplingNumNonzeros( SolverBase const * const solver,
                        DomainPartition & domain,
                        DofManager & dofManager,
                        arrayView1d< localIndex > const & rowLengths,
                        integer const resNumDof,
                        integer const wellNumDof,
                        string const & resElemDofName,
                        string const & wellElemDofName );

/**
 * @brief Validate the well perforations ensuring that each perforation is located in a reservoir region that is also
 * targetted by the solver
 * @param reservoirSolver the reservoir solver
 * @param wellSolver the well solver
 * @param domain the physical domain object
 */
bool validateWellPerforations( SolverBase const * const reservoirSolver,
                               WellSolverBase const * const wellSolver,
                               DomainPartition const & domain );

}

template< typename RESERVOIR_SOLVER, typename WELL_SOLVER >
class CoupledReservoirAndWellsBase : public CoupledSolver< RESERVOIR_SOLVER, WELL_SOLVER >
{
public:

  using Base = CoupledSolver< RESERVOIR_SOLVER, WELL_SOLVER >;
  using Base::m_solvers;
  using Base::m_names;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  enum class SolverType : integer
  {
    Reservoir = 0,
    Well = 1
  };

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "reservoirAndWells"; }

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  CoupledReservoirAndWellsBase ( const string & name,
                                 dataRepository::Group * const parent )
    : Base( name, parent ),
    m_isWellTransmissibilityComputed( false )
  {
    this->template getWrapper< string >( Base::viewKeyStruct::discretizationString() ).
      setInputFlag( dataRepository::InputFlags::FALSE );
  }

  /**
   * @brief default destructor
   */
  virtual ~CoupledReservoirAndWellsBase () override {}

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = true ) override
  {
    GEOS_MARK_FUNCTION;

    // call reservoir solver setup (needed in case of SinglePhasePoromechanicsConformingFractures)
    reservoirSolver()->setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

    dofManager.setDomain( domain );

    Base::setupDofs( domain, dofManager );
    dofManager.reorderByRank();

    // Set the sparsity pattern without reservoir-well coupling
    SparsityPattern< globalIndex > patternDiag;
    dofManager.setSparsityPattern( patternDiag );

    // Get the original row lengths (diagonal blocks only)
    array1d< localIndex > rowLengths( patternDiag.numRows() );
    for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
    {
      rowLengths[localRow] = patternDiag.numNonZeros( localRow );
    }

    // Add the number of nonzeros induced by coupling on perforations
    addCouplingNumNonzeros( domain, dofManager, rowLengths.toView() );

    // Create a new pattern with enough capacity for coupled matrix
    SparsityPattern< globalIndex > pattern;
    pattern.resizeFromRowCapacities< parallelHostPolicy >( patternDiag.numRows(), patternDiag.numColumns(), rowLengths.data() );

    // Copy the original nonzeros
    for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
    {
      globalIndex const * cols = patternDiag.getColumns( localRow ).dataIfContiguous();
      pattern.insertNonZeros( localRow, cols, cols + patternDiag.numNonZeros( localRow ) );
    }

    // Add the nonzeros from coupling
    addCouplingSparsityPattern( domain, dofManager, pattern.toView() );

    // Finally, steal the pattern into a CRS matrix
    localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
    localMatrix.setName( this->getName() + "/localMatrix" );

    rhs.setName( this->getName() + "/rhs" );
    rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

    solution.setName( this->getName() + "/solution" );
    solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );
  }

  /**@}*/

  /**
   * @brief accessor for the pointer to the reservoir solver
   * @return a pointer to the reservoir solver
   */
  RESERVOIR_SOLVER *
  reservoirSolver() const { return std::get< toUnderlying( SolverType::Reservoir ) >( m_solvers ); }

  /**
   * @brief accessor for the pointer to the well solver
   * @return a pointer to the well solver
   */
  WELL_SOLVER *
  wellSolver() const { return std::get< toUnderlying( SolverType::Well ) >( m_solvers ); }

  virtual void
  initializePostInitialConditionsPreSubGroups() override
  {
    Base::initializePostInitialConditionsPreSubGroups( );

    DomainPartition & domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );

    // Validate well perforations: Ensure that each perforation is in a region targeted by the solver
    if( !validateWellPerforations( domain ))
    {
      return;
    }
  }

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override
  {
    Base::implicitStepSetup( time_n, dt, domain );

    // we delay the computation of the transmissibility until the last minute
    // because we want to make sure that the permeability has been updated (in the flow solver)
    // this is necessary for some permeability models (like Karman-Kozeny) that do not use the imported permeability
    // ultimately, we may want to use this mechanism to update the well transmissibility at each time step (if needed)
    if( !m_isWellTransmissibilityComputed )
    {
      computeWellTransmissibility( domain );
      m_isWellTransmissibilityComputed = true;
    }
  }

  void
  assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const
  { reservoirSolver()->assembleFluxTerms( dt, domain, dofManager, localMatrix, localRhs );  }

  void
  assembleStabilizedFluxTerms( real64 const dt,
                               DomainPartition const & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) const
  { reservoirSolver()->assembleStabilizedFluxTerms( dt, domain, dofManager, localMatrix, localRhs );  }

  real64 updateFluidState( ElementSubRegionBase & subRegion ) const
  { return reservoirSolver()->updateFluidState( subRegion ); }
  void updatePorosityAndPermeability( CellElementSubRegion & subRegion ) const
  { reservoirSolver()->updatePorosityAndPermeability( subRegion ); }
  void updateSolidInternalEnergyModel( ObjectManagerBase & dataGroup ) const
  { reservoirSolver()->updateSolidInternalEnergyModel( dataGroup ); }

  integer & isThermal() { return reservoirSolver()->isThermal(); }

  void enableJumpStabilization()
  { reservoirSolver()->enableJumpStabilization(); }

  void enableFixedStressPoromechanicsUpdate()
  { reservoirSolver()->enableFixedStressPoromechanicsUpdate(); }

  void setKeepFlowVariablesConstantDuringInitStep( bool const keepFlowVariablesConstantDuringInitStep )
  { reservoirSolver()->setKeepFlowVariablesConstantDuringInitStep( keepFlowVariablesConstantDuringInitStep ); }

  virtual void saveSequentialIterationState( DomainPartition & domain ) override
  { reservoirSolver()->saveSequentialIterationState( domain ); }

protected:

  /**
   * @Brief loop over perforations and increase Jacobian matrix row lengths for reservoir and well elements accordingly
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param rowLengths the row-by-row length
   */
  void
  addCouplingNumNonzeros( DomainPartition & domain,
                          DofManager & dofManager,
                          arrayView1d< localIndex > const & rowLengths ) const
  {
    coupledReservoirAndWellsInternal::
      addCouplingNumNonzeros( this,
                              domain,
                              dofManager,
                              rowLengths,
                              wellSolver()->numDofPerResElement(),
                              wellSolver()->numDofPerWellElement(),
                              wellSolver()->resElementDofName(),
                              wellSolver()->wellElementDofName() );
  }

  /**
   * @Brief add the sparsity pattern induced by the perforations
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  virtual void
  addCouplingSparsityPattern( DomainPartition const & domain,
                              DofManager const & dofManager,
                              SparsityPatternView< globalIndex > const & pattern ) const = 0;

  /// Flag to determine whether the well transmissibility needs to be computed
  bool m_isWellTransmissibilityComputed;

private:

  /**
   * @brief Validate the well perforations ensuring that each perforation is located in a reservoir region that is also
   * targeted by the solver
   * @param domain the physical domain object
   */
  bool validateWellPerforations( DomainPartition const & domain ) const
  {
    return coupledReservoirAndWellsInternal::validateWellPerforations( reservoirSolver(), wellSolver(), domain );
  }

  /**
   * @brief Compute the transmissibility at all the perforations
   * @param domain the physical domain object
   */
  void computeWellTransmissibility( DomainPartition & domain ) const
  {
    this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                                 MeshLevel & meshLevel,
                                                                                 arrayView1d< string const > const & regionNames )
    {
      ElementRegionManager & elemManager = meshLevel.getElemManager();

      ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const elemCenter =
        elemManager.constructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

      // loop over the wells
      elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                  WellElementSubRegion & subRegion )
      {
        array1d< array1d< arrayView3d< real64 const > > > const permeability =
          elemManager.constructMaterialFieldAccessor< constitutive::PermeabilityBase,
                                                      fields::permeability::permeability >();

        PerforationData & perforationData = *subRegion.getPerforationData();
        WellControls const & wellControls = wellSolver()->getWellControls( subRegion );

        // compute the Peaceman index (if not read from XML)
        perforationData.computeWellTransmissibility( meshLevel, subRegion, permeability );

        // if the log level is 1, we output the value of the transmissibilities
        if( wellControls.getLogLevel() >= 2 )
        {
          arrayView2d< real64 const > const perfLocation =
            perforationData.getField< fields::perforation::location >();
          arrayView1d< real64 const > const perfTrans =
            perforationData.getField< fields::perforation::wellTransmissibility >();

          // get the element region, subregion, index
          arrayView1d< localIndex const > const resElemRegion =
            perforationData.getField< fields::perforation::reservoirElementRegion >();
          arrayView1d< localIndex const > const resElemSubRegion =
            perforationData.getField< fields::perforation::reservoirElementSubRegion >();
          arrayView1d< localIndex const > const resElemIndex =
            perforationData.getField< fields::perforation::reservoirElementIndex >();

          GEOS_UNUSED_VAR( perfLocation ); // unused if geos_error_if is nulld
          GEOS_UNUSED_VAR( perfTrans ); // unused if geos_error_if is nulld
          GEOS_UNUSED_VAR( resElemRegion ); // unused if geos_error_if is nulld
          GEOS_UNUSED_VAR( resElemSubRegion ); // unused if geos_error_if is nulld
          GEOS_UNUSED_VAR( resElemIndex ); // unused if geos_error_if is nulld

          forAll< serialPolicy >( perforationData.size(), [=] ( localIndex const iperf )
          {
            GEOS_UNUSED_VAR( iperf ); // unused if geos_error_if is nulld
            GEOS_LOG_RANK( GEOS_FMT( "Perforation at ({},{},{}); perforated element center: ({},{},{}); transmissibility: {} [{}]",
                                     perfLocation[iperf][0], perfLocation[iperf][1], perfLocation[iperf][2],
                                     elemCenter[resElemRegion[iperf]][resElemSubRegion[iperf]][resElemIndex[iperf]][0],
                                     elemCenter[resElemRegion[iperf]][resElemSubRegion[iperf]][resElemIndex[iperf]][1],
                                     elemCenter[resElemRegion[iperf]][resElemSubRegion[iperf]][resElemIndex[iperf]][2],
                                     perfTrans[iperf], getSymbol( units::Transmissibility ) ) );
          } );
        }
      } );
    } );
  }

};


} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDRESERVOIRANDWELLSBASE_HPP_ */
