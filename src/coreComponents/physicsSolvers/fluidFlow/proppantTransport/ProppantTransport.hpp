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
 * @file ProppantTransport.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORT_PROPPANTTRANSPORT_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORT_PROPPANTTRANSPORT_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "constitutive/fluid/singlefluid/ParticleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SlurryFluidBase.hpp"

namespace geos
{

namespace dataRepository
{
class Group;
}
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;

/**
 * @class ProppantTransport
 *
 * class to perform a proppant finite volume solve.
 */
class ProppantTransport : public FlowSolverBase
{
public:

  /// String used to form the solverName used to register single-physics solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "proppant"; }

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  ProppantTransport( const string & name,
                     Group * const parent );


  /// deleted default constructor
  ProppantTransport() = delete;

  /// deleted copy constructor
  ProppantTransport( ProppantTransport const & ) = delete;

  /// default move constructor
  ProppantTransport( ProppantTransport && ) = default;

  /// deleted assignment operator
  ProppantTransport & operator=( ProppantTransport const & ) = delete;

  /// deleted move operator
  ProppantTransport & operator=( ProppantTransport && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~ProppantTransport() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName() { return "ProppantTransport"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  void preStepUpdate( real64 const & time_n,
                      real64 const & dt,
                      DomainPartition & domain );

  void postStepUpdate( real64 const & time_n,
                       real64 const & dt,
                       DomainPartition & domain );

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
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time,
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
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  /**
   * @brief assembles the accumulation terms for all cells
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */

  void assembleAccumulationTerms( real64 const dt,
                                  DomainPartition const & domain,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs );

  /**
   * @brief assembles the flux terms for all cells
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */
  void assembleFluxTerms( real64 const dt,
                          DomainPartition const & domain,
                          DofManager const & dofManager,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs );

  /**@}*/

  void resizeFractureFields( MeshLevel & mesh, arrayView1d< string const > const & regionNames );

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr char const * proppantNamesString() { return "proppantNames"; }

    // these are used to store last converged time step values

    static constexpr char const * updateProppantPackingString() { return "updateProppantPacking"; }
    static constexpr char const * bridgingFactorString() { return "bridgingFactor"; }
    static constexpr char const * maxProppantConcentrationString() { return "maxProppantConcentration"; }
    static constexpr char const * proppantDiameterString() { return "proppantDiameter"; }
    static constexpr char const * proppantDensityString() { return "proppantDensity"; }
    static constexpr char const * criticalShieldsNumberString() { return "criticalShieldsNumber"; }
    static constexpr char const * frictionCoefficientString() { return "frictionCoefficient"; }
  };

  virtual void initializePostInitialConditionsPreSubGroups() override;

  void updateProppantMobility( ObjectManagerBase & dataGroup );

  /**
   * @brief Function to update proppant pack volume fraction
   */
  void updateProppantPackVolume( real64 const time_n, real64 const dt, DomainPartition & domain );

  virtual void updateState ( DomainPartition & domain ) override final { GEOS_UNUSED_VAR( domain ); };

  /**
   * @brief Function to update fluid and proppant properties
   * @param domain the domain
   */
  void updateState( ObjectManagerBase & dataGroup );

protected:

  virtual void postInputInitialization() override;

private:

  /**
   * @brief Function to update fluid properties
   * @param domain the domain
   */
  void updateFluidModel( ObjectManagerBase & dataGroup );

  void updateComponentDensity( ObjectManagerBase & dataGroup );

  void updateProppantModel( ObjectManagerBase & dataGroup );

  /**
   * @brief Function to update cell-based fluid flux
   */
  void updateCellBasedFlux( real64 const time_n,
                            DomainPartition & domain );

  void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;

  integer m_numComponents;

  integer m_updateProppantPacking;
  real64 m_bridgingFactor;
  real64 m_minAperture;
  real64 m_maxProppantConcentration;
  real64 m_proppantDiameter;
  real64 m_proppantDensity;
  real64 m_criticalShieldsNumber;
  real64 m_frictionCoefficient;
};


} /* namespace geos */

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORT_PROPPANTTRANSPORT_HPP_
