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
 * @file WellSolverBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASE_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASE_HPP_

#include "physicsSolvers/PhysicsSolverBase.hpp"

namespace geos
{

class DomainPartition;
class WellControls;
class WellElementSubRegion;
class TableFunction;

/**
 * @class WellSolverBase
 *
 * Base class for well solvers.
 * Provides some common features
 */
class WellSolverBase : public PhysicsSolverBase
{
public:

  /// String used to form the solverName used to register single-physics solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "well"; }

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  WellSolverBase( const string & name,
                  Group * const parent );

  /// default destructor
  virtual ~WellSolverBase() override;

  /// deleted default constructor
  WellSolverBase() = delete;

  /// deleted copy constructor
  WellSolverBase( WellSolverBase const & ) = delete;

  /// default move constructor
  WellSolverBase( WellSolverBase && ) = default;

  /// deleted assignment operator
  WellSolverBase & operator=( WellSolverBase const & ) = delete;

  /// deleted move operator
  WellSolverBase & operator=( WellSolverBase && ) = delete;

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// Expand catalog for schema generation
  virtual void expandObjectCatalogs() override;


  /**
   * @brief setter for the name of the flow solver (needed to use the flow kernels like UpdateFluid)
   * @param name the name of the flow solver
   */
  void setFlowSolverName( string const & name ) { m_flowSolverName = name; }

  /**
   * @brief getter for the name of the flow solver (used in UpdateState)
   * @return a string containing the name of the flow solver
   */
  string const & getFlowSolverName() const { return m_flowSolverName; }

  /**
   * @brief getter for the number of degrees of freedom per well element
   * @return the number of dofs
   */
  localIndex numDofPerWellElement() const { return m_numDofPerWellElement; }

  /**
   * @brief getter for the number of degrees of freedom per mesh element
   * @return the number of dofs
   */
  localIndex numDofPerResElement() const { return m_numDofPerResElement; }

  /**
   * @brief getter for iso/thermal switch
   * @return True if thermal
   */
  integer isThermal() const { return m_isThermal; }

  /**
   * @brief get the name of DOF defined on well elements
   * @return name of the DOF field used by derived solver type
   */
  virtual string wellElementDofName() const = 0;

  /**
   * @brief get the name of DOF defined on well elements
   * @return name of the DOF field used by derived solver type
   */
  virtual string resElementDofName() const = 0;

  /**
   * @brief const getter for the number of fluid components
   * @return the number of fluid components
   */
  virtual localIndex numFluidComponents() const = 0;

  /**
   * @brief Get the number of fluid phases
   * @return the number of phases
   */
  virtual localIndex numFluidPhases() const = 0;

  /**
   * @brief getter for the well controls associated to this well subRegion
   * @param subRegion the well subRegion whose controls are requested
   * @return a reference to the controls
   */
  WellControls & getWellControls( WellElementSubRegion const & subRegion );

  /**
   * @brief const getter for the well controls associated to this well subRegion
   * @param subRegion the well subRegion whose controls are requested
   * @return a reference to the const controls
   */
  WellControls const & getWellControls( WellElementSubRegion const & subRegion ) const;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain ) override;

  virtual void implicitStepComplete( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                     real64 const & GEOS_UNUSED_PARAM( dt ),
                                     DomainPartition & GEOS_UNUSED_PARAM( domain ) ) override {}

  virtual void applyBoundaryConditions( real64 const GEOS_UNUSED_PARAM( time_n ),
                                        real64 const GEOS_UNUSED_PARAM( dt ),
                                        DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                        DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                        CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                        arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) ) override {}


  /**@}*/

  /**
   * @brief function to assemble the linear system matrix and rhs
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

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
                                  arrayView1d< real64 > const & localRhs ) = 0;

  /**
   * @brief assembles the accumulation term for all the well elements
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleAccumulationTerms( real64 const & time_n,
                                          real64 const & dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs ) = 0;

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
                                          arrayView1d< real64 > const & localRhs ) = 0;

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param domain the domain containing the mesh and fields
   */
  virtual void updateState( DomainPartition & domain ) override;

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param subRegion the well subRegion containing the well elements and their associated fields
   */
  virtual real64 updateSubRegionState( WellElementSubRegion & subRegion ) = 0;

  /**
   * @brief Recompute the perforation rates for all the wells
   * @param domain the domain containing the mesh and fields
   */
  virtual void computePerforationRates( real64 const & time_n,
                                        real64 const & dt,
                                        DomainPartition & domain ) = 0;

  /**
   * @brief function to set the next time step size
   * @param[in] currentTime the current time
   * @param[in] lastDt the last time step size
   * @param[in] domain the domain object
   * @return the prescribed time step size
   */
  virtual real64 setNextDt( real64 const & currentTime,
                            real64 const & lastDt,
                            DomainPartition & domain ) override;

  /**
   * @brief Utility function to keep the well variables during a time step (used in poromechanics simulations)
   * @param[in] keepVariablesConstantDuringInitStep flag to tell the solver to freeze its primary variables during a time step
   * @detail This function is meant to be called by a specific task before/after the initialization step
   */
  void setKeepVariablesConstantDuringInitStep( bool const keepVariablesConstantDuringInitStep )
  { m_keepVariablesConstantDuringInitStep = keepVariablesConstantDuringInitStep; }

  struct viewKeyStruct : PhysicsSolverBase::viewKeyStruct
  {
    static constexpr char const * fluidNamesString() { return "fluidNames"; }
    static constexpr char const * isThermalString() { return "isThermal"; }
    static constexpr char const * writeCSVFlagString() { return "writeCSV"; }
    static constexpr char const * timeStepFromTablesFlagString() { return "timeStepFromTables"; }
  };

private:

  /**
   * @brief This function generates various discretization information for later use.
   * @param domain the domain parition
   */
  void precomputeData( DomainPartition & domain );

  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const override;


protected:

  virtual void postInputInitialization() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void initializePostSubGroups() override;

  /**
   * @brief Initialize all the primary and secondary variables in all the wells
   * @param domain the domain containing the well manager to access individual wells
   */
  virtual void initializeWells( DomainPartition & domain, real64 const & time_n, real64 const & dt ) = 0;

  /**
   * @brief Make sure that the well constraints are compatible
   * @param time_n the time at the beginning of the time step
   * @param dt the time step dt
   * @param subRegion the well subRegion
   */
  virtual void validateWellConstraints( real64 const & time_n,
                                        real64 const & dt,
                                        WellElementSubRegion const & subRegion ) = 0;

  virtual void printRates( real64 const & time_n,
                           real64 const & dt,
                           DomainPartition & domain ) = 0;

  void setNextDtFromTable( TableFunction const * table, real64 const currentTime, real64 & nextDt );

  /// name of the flow solver
  string m_flowSolverName;

  /// the max number of fluid phases
  integer m_numPhases;

  /// the number of fluid components
  integer m_numComponents;

  /// the number of Degrees of Freedom per well element
  integer m_numDofPerWellElement;

  /// the number of Degrees of Freedom per reservoir element
  integer m_numDofPerResElement;

  /// flag indicating whether thermal formulation is used
  integer m_isThermal;

  /// rates output
  integer m_writeCSV;
  string const m_ratesOutputDir;

  // flag to enable time step selection base on rates/bhp tables coordinates
  integer m_timeStepFromTables;

  /// flag to freeze the initial state during initialization in coupled problems
  integer m_keepVariablesConstantDuringInitStep;

  /// name of the fluid constitutive model used as a reference for component/phase description
  string m_referenceFluidModelName;
};

}

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASE_HPP_
