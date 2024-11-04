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
 * @file ReactiveCompositionalMultiphaseOBL.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_REACTIVECOMPOSITIONALMULTIPHASEOBL_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_REACTIVECOMPOSITIONALMULTIPHASEOBL_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "functions/MultivariableTableFunction.hpp"

namespace geos
{

/**
 * @class ReactiveCompositionalMultiphaseOBL
 * - A compositional multiphase thermal reactive solver based on OBL linearization approach (https://doi.org/10.1016/j.petrol.2017.08.009)
 * - Analogouos to Super Engine in DARTS reservoir simulator (https://darts.citg.tudelft.nl/)
 * - Takes into account diffusion, capillarity, kinetic and equilibrium reactions, precipitation/dissolution using only cell-centered
 * variables
 * - Supports only TPFA discretization
 * - Does not work with wells and aquifers (will require introduction of additional operator tables)
 * - Uses a single operator table for the whole reservoir (introduction of several tables will allow to support different fluid properties
 * in different reservoir regions)
 * - Since only static MultivariableTableFunction is currently implemented, all operator values have to be stored, which limits OBL
 * discretization:
 *    only 3-5 components with 32-64 points is now viable. Adative MultivariableTableFunction with sparse storage will allow for much more
 * refined OBL parametrizations
 * - Does not use any fluid model, and solid models are only needed to get initial porosity
 */
//START_SPHINX_INCLUDE_00
class ReactiveCompositionalMultiphaseOBL : public FlowSolverBase
{
//END_SPHINX_INCLUDE_00
public:


  /**
   * @brief Maximum supported number of fluid components (species)
   *
   * @note This puts an upper bound on memory use, allowing to optimize code better
   */
  static constexpr integer MAX_NUM_COMPONENTS = 5;

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  ReactiveCompositionalMultiphaseOBL( const string & name,
                                      Group * const parent );

  /// deleted default constructor
  ReactiveCompositionalMultiphaseOBL() = delete;

  /// deleted copy constructor
  ReactiveCompositionalMultiphaseOBL( ReactiveCompositionalMultiphaseOBL const & ) = delete;

  /// default move constructor
  ReactiveCompositionalMultiphaseOBL( ReactiveCompositionalMultiphaseOBL && ) = default;

  /// deleted assignment operator
  ReactiveCompositionalMultiphaseOBL & operator=( ReactiveCompositionalMultiphaseOBL const & ) = delete;

  /// deleted move operator
  ReactiveCompositionalMultiphaseOBL & operator=( ReactiveCompositionalMultiphaseOBL && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~ReactiveCompositionalMultiphaseOBL() override = default;

//START_SPHINX_INCLUDE_01
  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string catalogName() { return "ReactiveCompositionalMultiphaseOBL"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Getter for the fluid component names
   * @return an array storing the component names
   */
  arrayView1d< string const > componentNames() const { return m_componentNames; }

  virtual void registerDataOnMesh( Group & meshBodies ) override;
//END_SPHINX_INCLUDE_01

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
  assembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

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

  virtual void updateState( DomainPartition & domain ) override final;


  /**
   * @brief Recompute operator values and derivatives from primary variables
   * @param dataGroup the group storing the required fields
   */
  void updateOBLOperators( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Get the number of fluid components (species)
   * @return the number of components
   */
  localIndex numFluidComponents() const { return m_numComponents; }

  /**
   * @brief Get the number of fluid phases
   * @return the number of phases
   */
  localIndex numFluidPhases() const { return m_numPhases; }

  /**
   * @brief Get the number of OBL operators
   * @return the number of OBL operators
   */
  localIndex numOBLOperators() const { return m_numOBLOperators; }

  /**
   * @brief assembles the accumulation term for all cells
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  void assembleAccumulationTerms( real64 const dt,
                                  DomainPartition & domain,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief assembles the flux term for all cells
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  virtual void
  assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const;

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr char const * elemDofFieldString() { return "compositionalMolarVariables"; }

    // inputs

    static constexpr char const * numPhasesString() { return "numPhases"; }

    static constexpr char const * numComponentsString() { return "numComponents"; }

    static constexpr char const * enableEnergyBalanceString() { return "enableEnergyBalance"; }

    static constexpr char const * componentNamesString() { return "componentNames"; }

    static constexpr char const * phaseNamesString() { return "phaseNames"; }

    static constexpr char const * transMultExpString() { return "transMultExp"; }

    static constexpr char const * maxCompFracChangeString() { return "maxCompFractionChange"; }

    static constexpr char const * allowLocalOBLChoppingString() { return "allowLocalOBLChopping"; }

    static constexpr char const * useDARTSL2NormString() { return "useDARTSL2Norm"; }
  };


  /**
   * @brief Function to perform checks before applying Dirichlet type BC's
   * @param domain the domain
   * @param time the current time
   */
  bool validateDirichletBC( DomainPartition & domain,
                            real64 const time ) const;

  /**
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void applyDirichletBC( real64 const time,
                         real64 const dt,
                         DofManager const & dofManager,
                         DomainPartition & domain,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Apply source flux boundary conditions to the system
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void applySourceFluxBC( real64 const time,
                          real64 const dt,
                          DofManager const & dofManager,
                          DomainPartition & domain,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs ) const;


  /**
   * @brief Sets all the primary variables, which are outside of defined OBL parametrization limits, to be inside them
   * Helps to evaluate correct operator values/derivatives when a disturbed solution is obtained from a previous Newton iteration
   * Does not help in case if the model is defined inconsistently:
   *    for instance, if minimum pressure value in OBL table is higher, than BHP used by a producer/sink
   *    In that case, the solution will be chopped all the time and Newton will never converge
   * @param domain the physical domain object
   */
  void chopPrimaryVariables( DomainPartition & domain );

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void initializePreSubGroups() override;


private:

  virtual void postInputInitialization() override;

  /// the max number of fluid phases
  integer m_numPhases;

  /// the number of fluid components
  integer m_numComponents;

  /// list of component names
  array1d< string > m_componentNames;

  /// list of phase names names
  array1d< string > m_phaseNames;

  /// the number of OBL operators
  integer m_numOBLOperators;

  /// flag indicating whether energy balance will be enabled or not
  integer m_enableEnergyBalance;

  /// maximum (absolute) change in a component fraction between two Newton iterations
  real64 m_maxCompFracChange;

  /// exponent of dynamic transmissibility multiplier
  real64 m_transMultExp;

  /// minimum value of the scaling factor obtained by enforcing maxCompFracChange
  real64 m_minScalingFactor;

  /// flag indicating whether local (cell-wise) chopping of negative compositions is allowed
  integer m_allowOBLChopping;

  /// flag indicating whether DARTS L2 norm is used for Newton convergence criterion
  integer m_useDARTSL2Norm;
};


} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_REACTIVECOMPOSITIONALMULTIPHASEOBL_HPP_
