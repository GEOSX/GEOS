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
 * @file SinglePhaseHybridFVM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVM_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVM_HPP_

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geos
{


/**
 * @class SinglePhaseHybridFVM
 *
 * class to assemble the single-phase flow equations based
 * on a mixed hybrid formulation known as the mimetic method
 */
class SinglePhaseHybridFVM : public SinglePhaseBase
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseHybridFVM( const string & name,
                        Group * const parent );


  /// deleted default constructor
  SinglePhaseHybridFVM() = delete;

  /// deleted copy constructor
  SinglePhaseHybridFVM( SinglePhaseHybridFVM const & ) = delete;

  /// default move constructor
  SinglePhaseHybridFVM( SinglePhaseHybridFVM && ) = default;

  /// deleted assignment operator
  SinglePhaseHybridFVM & operator=( SinglePhaseHybridFVM const & ) = delete;

  /// deleted move operator
  SinglePhaseHybridFVM & operator=( SinglePhaseHybridFVM && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseHybridFVM() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  { return "SinglePhaseHybridFVM"; }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

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

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) override;

  virtual void
  assembleStabilizedFluxTerms( real64 const dt,
                               DomainPartition const & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void
  assembleEDFMFluxTerms( real64 const time_n,
                         real64 const dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         string const & jumpDofKey ) override final;

  virtual void
  assembleHydrofracFluxTerms( real64 const time_n,
                              real64 const dt,
                              DomainPartition const & domain,
                              DofManager const & dofManager,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs,
                              CRSMatrixView< real64, localIndex const > const & dR_dAper ) override final;

  virtual void
  applyAquiferBC( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) const override final;

  virtual void
  saveAquiferConvergedState( real64 const & time,
                             real64 const & dt,
                             DomainPartition & domain ) override;

  virtual void updatePressureGradient( DomainPartition & domain ) override final;

  /**
   * @brief Function to perform the application of Dirichlet BCs on faces
   * @param[in] time_n current time
   * @param[in] dt time step
   * @param[in] faceSet degree-of-freedom manager associated with the linear system
   * @param[in] domain the domain
   * @param[inout] localMatrix the system matrix
   * @param[inout] localRhs the system right-hand side vector
   */
  void
  applyFaceDirichletBC( real64 const time_n,
                        real64 const dt,
                        DofManager const & faceSet,
                        DomainPartition & domain,
                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                        arrayView1d< real64 > const & localRhs );

  /**@}*/

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

private:

  /// relative tolerance (redundant with FluxApproximationBase)
  real64 m_areaRelTol;

  /// region filter used in flux assembly
  SortedArray< localIndex > m_regionFilter;

};

} /* namespace geos */

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVM_HPP_
