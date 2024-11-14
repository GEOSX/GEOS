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
 * @file ImmiscibleMultiphaseFlow.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOW_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOW_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/ImmiscibleMultiphaseKernels.hpp"

namespace geos
{

//START_SPHINX_INCLUDE_00
/**
 * @class ImmiscibleMultiphaseFlow
 *
 * An Immiscible multiphase flow solver
 */
class ImmiscibleMultiphaseFlow : public FlowSolverBase
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  ImmiscibleMultiphaseFlow( const string & name,
                            Group * const parent );

  /// deleted default constructor
  ImmiscibleMultiphaseFlow() = delete;

  /// deleted copy constructor
  ImmiscibleMultiphaseFlow( ImmiscibleMultiphaseFlow const & ) = delete;

  /// default move constructor
  ImmiscibleMultiphaseFlow( ImmiscibleMultiphaseFlow && ) = default;

  /// deleted assignment operator
  ImmiscibleMultiphaseFlow & operator=( ImmiscibleMultiphaseFlow const & ) = delete;

  /// deleted move operator
  ImmiscibleMultiphaseFlow & operator=( ImmiscibleMultiphaseFlow && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~ImmiscibleMultiphaseFlow() override = default;
  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string catalogName() { return "ImmiscibleMultiphaseFlow"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( Group & meshBodies ) override final;

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
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

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

  void updateFluidState( ElementSubRegionBase & subRegion ) const;

  void updatePhaseVolumeFraction( ElementSubRegionBase & subRegion ) const;

  virtual void saveConvergedState( ElementSubRegionBase & subRegion ) const override final;

  virtual void updateState( DomainPartition & domain ) override final;

  /**
   * @brief Getter for the number of fluid phases
   * @return the number of phases
   */
  integer numFluidPhases() const { return m_numPhases; }

  /**
   * @brief assembles the accumulation and volume balance terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  void assembleAccumulationTerm( DomainPartition & domain,
                                 DofManager const & dofManager,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const;
  /**
   * @brief Initialize all variables from initial conditions
   * @param domain the domain containing the mesh and fields
   *
   * Initialize all variables from initial conditions. This calculating primary variable values
   * from prescribed intermediate values (i.e. global densities from global fractions)
   * and any applicable hydrostatic equilibration of the domain
   */
  void initializeFluidState( MeshLevel & mesh, DomainPartition & domain, arrayView1d< string const > const & regionNames );

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

  void applySourceFluxBC( real64 const time,
                          real64 const dt,
                          DofManager const & dofManager,
                          DomainPartition & domain,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs ) const;
  /**
   * @brief function to set the next time step size
   * @param[in] currentDt the current time step size
   * @param[in] domain the domain object
   * @return the prescribed time step size
   */
  real64 setNextDt( real64 const & currentDt,
                    DomainPartition & domain ) override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  /**
   * @brief Function to update fluid mass
   * @param subRegion subregion that contains the fields
   */
  void updatePhaseMass( ElementSubRegionBase & subRegion ) const;

  struct viewKeyStruct : public FlowSolverBase::viewKeyStruct
  {
    static constexpr char const * capPressureNamesString() { return "capPressureNames"; }
    static constexpr char const * relPermNamesString() { return "relPermNames"; }
    static constexpr char const * elemDofFieldString() { return "elemDofField"; }
    static constexpr char const * useTotalMassEquationString() { return "useTotalMassEquation"; }
  };


private:

  virtual void postInputInitialization() override;

  virtual void initializePreSubGroups() override;

  /**
   * @brief Update all relevant fluid models using current values of pressure and phase volume fraction
   * @param dataGroup the group storing the required fields
   */
  void updateFluidModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Update all relevant relperm models using current values of phase volume fraction
   * @param dataGroup the group storing the required fields
   */
  void updateRelPermModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Update all relevant capillary pressure models using current values of phase volume fraction
   * @param dataGroup the group storing the required fields
   */
  void updateCapPressureModel( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Recompute phase mobility from constitutive and primary variables
   * @param dataGroup the group storing the required field
   */
  void updatePhaseMobility( ObjectManagerBase & dataGroup ) const;

  /**
   * @brief Utility function that encapsulates the call to FieldSpecificationBase::applyFieldValue in BC application
   * @param[in] time_n the time at the beginning of the step
   * @param[in] dt the time step
   * @param[in] mesh the mesh level object
   * @param[in] logMessage the log message issued by the solver if the bc is called
   * @param[in] fieldKey the key of the field specified in the xml file
   * @param[in] boundaryFieldKey the key of the boundary field
   */
  template< typename OBJECT_TYPE >
  void applyFieldValue( real64 const & time_n,
                        real64 const & dt,
                        MeshLevel & mesh,
                        char const logMessage[],
                        string const fieldKey,
                        string const boundaryFieldKey ) const;

  /// the max number of fluid phases
  integer m_numPhases;

  /// flag to determine whether or not to apply capillary pressure
  integer m_hasCapPressure;

  /// flag to determine whether or not to use total velocity formulation
  integer m_useTotalMassEquation;



private:

  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;

};

template< typename OBJECT_TYPE >
void ImmiscibleMultiphaseFlow::applyFieldValue( real64 const & time_n,
                                                real64 const & dt,
                                                MeshLevel & mesh,
                                                char const logMessage[],
                                                string const fieldKey,
                                                string const boundaryFieldKey ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  fsManager.apply< OBJECT_TYPE >( time_n + dt,
                                  mesh,
                                  fieldKey,
                                  [&]( FieldSpecificationBase const & fs,
                                       string const & setName,
                                       SortedArrayView< localIndex const > const & lset,
                                       OBJECT_TYPE & targetGroup,
                                       string const & )
  {
    if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
    {
      globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( lset.size() );
      GEOS_LOG_RANK_0( GEOS_FMT( logMessage,
                                 getName(), time_n+dt, fs.getCatalogName(), fs.getName(),
                                 setName, targetGroup.getName(), fs.getScale(), numTargetElems ) );
    }

    // Specify the bc value of the field
    fs.applyFieldValue< FieldSpecificationEqual,
                        parallelDevicePolicy<> >( lset,
                                                  time_n + dt,
                                                  targetGroup,
                                                  boundaryFieldKey );
  } );
}


} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOW_HPP_
