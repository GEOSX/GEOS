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
 * @file SinglePhaseFVM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVM_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVM_HPP_

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseProppantBase.hpp"

namespace geos
{


/**
 * @class SinglePhaseFVM
 *
 * class to perform a single phase finite volume solve
 * using only cell-centered variables
 * works with both TPFA and MPFA
 */
template< typename BASE = SinglePhaseBase >
class SinglePhaseFVM : public BASE
{
public:


  // Aliasing public/protected members/methods of Group so we don't
  // have to use this->member etc.
  using BASE::getLogLevel;

  // Aliasing public/protected members/methods of SolverBase so we don't
  // have to use this->member etc.
  using BASE::forDiscretizationOnMeshTargets;
  using BASE::m_cflFactor;
  using BASE::m_maxStableDt;
  using BASE::m_nextDt;
  using BASE::m_discretizationName;
  using BASE::m_dofManager;
  using BASE::m_matrix;
  using BASE::m_rhs;
  using BASE::m_solution;
  using BASE::m_localMatrix;
  using BASE::m_linearSolverParameters;
  using BASE::m_nonlinearSolverParameters;

  // Aliasing public/protected members/methods of FlowSolverBase so we don't
  // have to use this->member etc.
  using BASE::m_numDofPerCell;
  using BASE::m_isThermal;

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseFVM( const string & name,
                  dataRepository::Group * const parent );


  /// deleted default constructor
  SinglePhaseFVM() = delete;

  /// deleted copy constructor
  SinglePhaseFVM( SinglePhaseFVM const & ) = delete;

  /// default move constructor
  SinglePhaseFVM( SinglePhaseFVM && ) = default;

  /// deleted assignment operator
  SinglePhaseFVM & operator=( SinglePhaseFVM const & ) = delete;

  /// deleted move operator
  SinglePhaseFVM & operator=( SinglePhaseFVM && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseFVM() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< BASE, SinglePhaseBase > )
    {
      return "SinglePhaseFVM";
    }
    else if constexpr ( std::is_same_v< BASE, SinglePhaseProppantBase > )
    {
      return "SinglePhaseProppantFVM";
    }
    else
    {
      return BASE::catalogName();
    }
  }

  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /**
   * @defgroup Solver Interface Function
   *
   * This function provides the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = true ) override;

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

  /**@}*/

  virtual void
  applyAquiferBC( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) const override;

  virtual void initializePreSubGroups() override;

private:

  /**
   * @brief Function to perform the application of Dirichlet BCs on faces
   * @param time_n current time
   * @param dt time step
   * @param faceSet degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void applyFaceDirichletBC( real64 const time_n,
                             real64 const dt,
                             DofManager const & faceSet,
                             DomainPartition & domain,
                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                             arrayView1d< real64 > const & localRhs );

  // no data needed here, see SinglePhaseBase

};

} /* namespace geos */

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVM_HPP_
