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
 * @file CompositionalMultiphaseHybridFVM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"

namespace geos
{

/**
 * @class CompositionalMultiphaseHybridFVM
 *
 * A compositional multiphase solver
 * using only cell-centered and face-centered variables
 */
class CompositionalMultiphaseHybridFVM : public CompositionalMultiphaseBase
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  CompositionalMultiphaseHybridFVM( const string & name,
                                    Group * const parent );

  /// deleted default constructor
  CompositionalMultiphaseHybridFVM() = delete;

  /// deleted copy constructor
  CompositionalMultiphaseHybridFVM( CompositionalMultiphaseHybridFVM const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseHybridFVM( CompositionalMultiphaseHybridFVM && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseHybridFVM & operator=( CompositionalMultiphaseHybridFVM const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseHybridFVM & operator=( CompositionalMultiphaseHybridFVM && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseHybridFVM() override = default;

  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string catalogName() { return "CompositionalMultiphaseHybridFVM"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( Group & MeshBodies ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( real64 const & time_n,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual real64
  scalingForSystemSolution( DomainPartition & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override;

  virtual bool
  checkSystemSolution( DomainPartition & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const override;

  virtual void
  assembleStabilizedFluxTerms( real64 const dt,
                               DomainPartition const & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) const override;

  virtual void
  updatePhaseMobility( ObjectManagerBase & dataGroup ) const override;

  virtual void
  applyAquiferBC( real64 const time,
                  real64 const dt,
                  DofManager const & dofManager,
                  DomainPartition & domain,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) const override;

  virtual void
  saveAquiferConvergedState( real64 const & time,
                             real64 const & dt,
                             DomainPartition & domain ) override;

  /**@}*/

  struct viewKeyStruct : CompositionalMultiphaseBase::viewKeyStruct
  {
    static constexpr char const * faceDofFieldString() { return "faceCenteredVariables"; }
  };

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void initializePreSubGroups() override;

protected:

  /// precompute the minGravityCoefficient for the buoyancy term
  void precomputeData( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) override;

private:

  /// tolerance used in the  computation of the transmissibility matrix
  real64 m_lengthTolerance;

  /// region filter used in flux assembly
  SortedArray< localIndex > m_regionFilter;

};

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_
