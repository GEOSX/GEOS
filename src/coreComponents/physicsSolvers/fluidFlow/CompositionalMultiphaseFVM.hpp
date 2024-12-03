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
 * @file CompositionalMultiphaseFVM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVM_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVM_HPP_

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"

namespace geos
{

/**
 * @class CompositionalMultiphaseFVM
 *
 * A compositional multiphase solver
 * using only cell-centered variables
 * works with both TPFA and MPFA
 */
//START_SPHINX_INCLUDE_00
class CompositionalMultiphaseFVM : public CompositionalMultiphaseBase
{
//END_SPHINX_INCLUDE_00
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  CompositionalMultiphaseFVM( const string & name,
                              Group * const parent );

  /// deleted default constructor
  CompositionalMultiphaseFVM() = delete;

  /// deleted copy constructor
  CompositionalMultiphaseFVM( CompositionalMultiphaseFVM const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseFVM( CompositionalMultiphaseFVM && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseFVM & operator=( CompositionalMultiphaseFVM const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseFVM & operator=( CompositionalMultiphaseFVM && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseFVM() override = default;

//START_SPHINX_INCLUDE_01
  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string catalogName() { return "CompositionalMultiphaseFVM"; }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }
//END_SPHINX_INCLUDE_01

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

  /**@}*/

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

  struct viewKeyStruct : CompositionalMultiphaseBase::viewKeyStruct
  {
    // DBC parameters
    static constexpr char const * useDBCString()                  { return "useDBC"; }
    static constexpr char const * omegaDBCString()                { return "omegaDBC"; }
    static constexpr char const * continuationDBCString()         { return "continuationDBC"; }

    static constexpr char const * miscibleDBCString()             { return "miscibleDBC"; }
    static constexpr char const * kappaminDBCString()             { return "kappaminDBC"; }
    static constexpr char const * contMultiplierDBCString()       { return "contMultiplierDBC"; }

    // nonlinear solver parameters
    static constexpr char const * scalingTypeString()               { return "scalingType"; }
  };

  /**
   * @brief Solution scaling type
   */
  enum class ScalingType : integer
  {
    Global,         ///< Scale the Newton update with a unique scaling factor
    Local            ///< Scale the Newton update locally (modifies the Newton direction)
  };

  /**
   * @brief Storage for value and element location, used to determine global max + location
   */
  template< typename VALUE_TYPE, typename INDEX_TYPE >
  struct valueAndLocation
  {
    valueAndLocation(){}
    valueAndLocation( VALUE_TYPE val, INDEX_TYPE loc ): value( val ), location( loc ){}
    VALUE_TYPE value;
    INDEX_TYPE location;
  };
  typedef valueAndLocation< real64, globalIndex > valueAndLocationType;

protected:

  virtual void postInputInitialization() override;

  virtual void
  initializePreSubGroups() override;

  struct DBCParameters
  {
    /// Flag to enable Dissipation Based Continuation Method
    integer useDBC;
    /// Factor by which the DBC flux is multiplied
    real64 omega;
    /// Factor by which the DBC flux is diminished every Newton
    real64 kappa;
    /// Flag to enable continuation for DBC Method
    integer continuation;
    /// Flag to enable DBC formulation
    integer miscible;
    /// Factor that controls how much dissipation is kept in the system when continuation is used
    real64 kappamin;
    /// Factor by which continuation parameter is changed every newton when DBC is used
    real64 contMultiplier;
  } m_dbcParams;

  /// Solution scaling type
  ScalingType m_scalingType;

private:

  /**
   * @brief Utility function to validate the consistency of face Dirichlet BC input
   * @param[in] domain the domain partition
   * @param[in] time the time at the end of the time step (time_n + dt)
   */
  bool validateFaceDirichletBC( DomainPartition & domain,
                                real64 const time ) const;

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

  // no data needed here, see CompositionalMultiphaseBase

};

ENUM_STRINGS( CompositionalMultiphaseFVM::ScalingType,
              "Global",
              "Local" );

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVM_HPP_
