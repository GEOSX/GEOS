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
 * @file SinglePhaseWell.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELL_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELL_HPP_

#include "WellSolverBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geos
{

namespace dataRepository
{
class Group;
}

namespace constitutive
{
class SingleFluidBase;
}
class WellElementSubRegion;

/**
 * @class SinglePhaseWell
 *
 * A single-phase well solver
 */
class SinglePhaseWell : public WellSolverBase
{
public:


  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseWell( const string & name,
                   Group * const parent );

  /// deleted default constructor
  SinglePhaseWell() = delete;

  /// deleted copy constructor
  SinglePhaseWell( SinglePhaseWell const & ) = delete;

  /// default move constructor
  SinglePhaseWell( SinglePhaseWell && ) = default;

  /// deleted assignment operator
  SinglePhaseWell & operator=( SinglePhaseWell const & ) = delete;

  /// deleted move operator
  SinglePhaseWell & operator=( SinglePhaseWell && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseWell() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName() { return "SinglePhaseWell"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual real64
  calculateResidualNorm( real64 const & time_n,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

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
  implicitStepSetup( real64 const & time,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  /**@}*/

  virtual string wellElementDofName() const override { return viewKeyStruct::dofFieldString(); }

  virtual string resElementDofName() const override;

  virtual localIndex numFluidComponents() const override { return 1; }

  virtual localIndex numFluidPhases() const override { return 1; }

  /**
   * @brief Recompute the volumetric rate that are used in the well constraints
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  virtual void updateVolRateForConstraint( WellElementSubRegion & subRegion );

  /**
   * @brief Recompute the BHP pressure that is used in the well constraints
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  virtual void updateBHPForConstraint( WellElementSubRegion & subRegion );

  /**
   * @brief Update fluid constitutive model state
   * @param subRegion the well subRegion containing the well elements and their associated fields
   */
  virtual void updateFluidModel( WellElementSubRegion & subRegion ) const;

  /**
   * @brief Recompute the perforation rates for all the wells
   * @param domain the domain containing the mesh and fields
   */
  virtual void computePerforationRates( real64 const & time_n,
                                        real64 const & dt, DomainPartition & domain ) override;

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models) on the well
   * @param subRegion the well subRegion containing the well elements and their associated fields
   */
  virtual real64 updateSubRegionState( WellElementSubRegion & subRegion ) override;

  /**
   * @brief assembles the flux terms for all connections between well elements
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleFluxTerms( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs ) override;

  /**
   * @brief assembles the accumulation term for all the well elements
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleAccumulationTerms( real64 const & time_n,
                                          real64 const & dt, DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs ) override;

  /**
   * @brief assembles the volume balance terms for all well elements
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void assembleVolumeBalanceTerms( DomainPartition const & domain,
                                   DofManager const & dofManager,
                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                   arrayView1d< real64 > const & localRhs );

  /**
   * @brief assembles the pressure relations at all connections between well elements except at the well head
   * @param time_n time at the beginning of the time step
   * @param dt the time step size
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assemblePressureRelations( real64 const & time_n,
                                          real64 const & dt,
                                          DomainPartition const & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs ) override;

  /*
   * @brief apply a special treatment to the wells that are shut
   * @param time_n the time at the previous converged time step
   * @param dt the time step size
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void shutDownWell( real64 const time_n,
                             real64 const dt,
                             DomainPartition const & domain,
                             DofManager const & dofManager,
                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                             arrayView1d< real64 > const & localRhs ) override;


  struct viewKeyStruct : WellSolverBase::viewKeyStruct
  {
    static constexpr char const * dofFieldString() { return "singlePhaseWellVars"; }

    // control data (not registered on the mesh)
    static constexpr char const * currentBHPString() { return "currentBHP"; }
    static constexpr char const * dCurrentBHP_dPresString() { return "dCurrentBHP_dPres"; }

    static constexpr char const * currentVolRateString() { return "currentVolumetricRate"; }
    static constexpr char const * dCurrentVolRate_dPresString() { return "dCurrentVolumetricRate_dPres"; }
    static constexpr char const * dCurrentVolRate_dRateString() { return "dCurrentVolumetricRate_dRate"; }

  };

protected:

  void printRates( real64 const & time_n,
                   real64 const & dt,
                   DomainPartition & domain ) override;

  /// flag if negative pressure is allowed
  integer m_allowNegativePressure;

private:

  /**
   * @brief Initialize all the primary and secondary variables in all the wells
   * @param domain the domain containing the well manager to access individual wells
   */
  void initializeWells( DomainPartition & domain, real64 const & time_n, real64 const & dt ) override;

  /**
   * @brief Make sure that the well constraints are compatible
   * @param time_n the time at the beginning of the time step
   * @param dt the time step dt
   * @param subRegion the well subRegion
   */
  virtual void validateWellConstraints( real64 const & time_n,
                                        real64 const & dt,
                                        WellElementSubRegion const & subRegion ) override;

};

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELL_HPP_
