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
 * @file CompositionalMultiphaseWell.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_

#include "physicsSolvers/fluidFlow/wells/WellSolverBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"

namespace geos
{

namespace constitutive
{
class ConstitutiveManager;
class MultiFluidBase;
}

/**
 * @class CompositionalMultiphaseWell
 *
 * A compositional multiphase well solver
 */
class CompositionalMultiphaseWell : public WellSolverBase
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  CompositionalMultiphaseWell( const string & name,
                               Group * const parent );

  /// deleted default constructor
  CompositionalMultiphaseWell() = delete;

  /// deleted copy constructor
  CompositionalMultiphaseWell( CompositionalMultiphaseWell const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseWell( CompositionalMultiphaseWell && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseWell & operator=( CompositionalMultiphaseWell const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseWell & operator=( CompositionalMultiphaseWell && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseWell() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName() { return "CompositionalMultiphaseWell"; }
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
  implicitStepSetup( real64 const & time,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  /**@}*/

  /**
   * @brief Recompute global component fractions from primary variables (component densities)
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void updateGlobalComponentFraction( WellElementSubRegion & subRegion ) const;

  /**
   * @brief Recompute the volumetric rates that are used in the well constraints
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */
  void updateVolRatesForConstraint( WellElementSubRegion & subRegion );

  /**
   * @brief Recompute the current BHP pressure
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */
  void updateBHPForConstraint( WellElementSubRegion & subRegion );

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */
  void updateFluidModel( WellElementSubRegion & subRegion );

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param subRegion the well subregion containing all the primary and dependent fields
   * @param targetIndex the targetIndex of the subRegion
   */
  void updatePhaseVolumeFraction( WellElementSubRegion & subRegion ) const;

  /**
   * @brief Recompute total mass densities from mass density and phase volume fractions
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void updateTotalMassDensity( WellElementSubRegion & subRegion ) const;

  /**
   * @brief Recompute the perforation rates for all the wells
   * @param domain the domain containing the mesh and fields
   */
  virtual void computePerforationRates( DomainPartition & domain ) override;

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  virtual void updateSubRegionState( WellElementSubRegion & subRegion ) override;

  virtual string wellElementDofName() const override { return viewKeyStruct::dofFieldString(); }

  virtual string resElementDofName() const override { return CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString(); }

  virtual localIndex numFluidComponents() const override { return m_numComponents; }

  virtual localIndex numFluidPhases() const override { return m_numPhases; }

  /**
   * @brief assembles the flux terms for all connections between well elements
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void assembleFluxTerms( real64 const dt,
                                  DomainPartition const & domain,
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
  virtual void assembleAccumulationTerms( DomainPartition const & domain,
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
  virtual void assembleVolumeBalanceTerms( DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs ) override;

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

  /**
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

  /**
   * @brief Sets all the negative component densities (if any) to zero.
   * @param domain the physical domain object
   */
  void chopNegativeDensities( DomainPartition & domain );

  arrayView1d< string const > relPermModelNames() const { return m_relPermModelNames; }

  struct viewKeyStruct : WellSolverBase::viewKeyStruct
  {
    static constexpr char const * dofFieldString() { return "compositionalWellVars"; }

    // inputs

    static constexpr char const * useMassFlagString() { return CompositionalMultiphaseBase::viewKeyStruct::useMassFlagString(); }

    static constexpr char const * useTotalMassEquationString() { return CompositionalMultiphaseBase::viewKeyStruct::useTotalMassEquationString(); }

    static constexpr char const * relPermNamesString() { return CompositionalMultiphaseBase::viewKeyStruct::relPermNamesString(); }

    static constexpr char const * maxCompFracChangeString() { return CompositionalMultiphaseBase::viewKeyStruct::maxCompFracChangeString(); }

    static constexpr char const * maxRelativePresChangeString() { return "maxRelativePressureChange"; }

    static constexpr char const * maxAbsolutePresChangeString() { return "maxAbsolutePressureChange"; }

    static constexpr char const * maxRelativeCompDensChangeString() { return "maxRelativeCompDensChange"; }

    static constexpr char const * allowLocalCompDensChoppingString() { return CompositionalMultiphaseBase::viewKeyStruct::allowLocalCompDensChoppingString(); }

    // control data (not registered on the mesh)

    static constexpr char const * massDensityString() { return "massDensity";}
    static constexpr char const * currentBHPString() { return "currentBHP"; }

    static constexpr char const * dCurrentBHP_dPresString() { return "dCurrentBHP_dPres"; }

    static constexpr char const * dCurrentBHP_dCompDensString() { return "dCurrentBHP_dCompDens"; }

    static constexpr char const * currentPhaseVolRateString() { return "currentPhaseVolumetricRate"; }

    static constexpr char const * dCurrentPhaseVolRate_dPresString() { return "dCurrentPhaseVolumetricRate_dPres"; }

    static constexpr char const * dCurrentPhaseVolRate_dCompDensString() { return "dCurrentPhaseVolumetricRate_dCompDens"; }

    static constexpr char const * dCurrentPhaseVolRate_dRateString() { return "dCurrentPhaseVolumetricRate_dRate"; }

    static constexpr char const * currentTotalVolRateString() { return "currentTotalVolumetricRate"; }

    static constexpr char const * dCurrentTotalVolRate_dPresString() { return "dCurrentTotalVolumetricRate_dPres"; }

    static constexpr char const * dCurrentTotalVolRate_dCompDensString() { return "dCurrentTotalVolumetricRate_dCompDens"; }

    static constexpr char const * dCurrentTotalVolRate_dRateString() { return "dCurrentTotalVolumetricRate_dRate"; }

  } viewKeysCompMultiphaseWell;

protected:

  virtual void postInputInitialization() override;

  virtual void initializePostSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  /*
   * @brief Utility function that checks the consistency of the constitutive models
   * @param[in] domain the domain partition
   * @detail
   * This function will produce an error if one of the well constitutive models
   * (fluid, relperm) is incompatible with the corresponding models in reservoir
   * regions connected to that particular well.
   */
  void validateConstitutiveModels( DomainPartition const & domain ) const;

  /**
   * @brief Checks if the WellControls parameters are within the fluid tables ranges
   * @param fluid the fluid to check
   */
  void validateWellControlsForFluid( WellControls const & wellControls,
                                     constitutive::MultiFluidBase const & fluid ) const;

  /**
   * @brief Checks injection streams for validity (compositions sum to one)
   * @param subRegion the well subRegion
   */
  void validateInjectionStreams( WellElementSubRegion const & subRegion ) const;

  /**
   * @brief Make sure that the well constraints are compatible
   * @param time_n the time at the beginning of the time step
   * @param dt the time step dt
   * @param subRegion the well subRegion
   */
  virtual void validateWellConstraints( real64 const & time_n,
                                        real64 const & dt,
                                        WellElementSubRegion const & subRegion ) override;

  void printRates( real64 const & time_n,
                   real64 const & dt,
                   DomainPartition & domain ) override;

private:

  /**
   * @brief Initialize all the primary and secondary variables in all the wells
   * @param domain the domain containing the well manager to access individual wells
   */
  void initializeWells( DomainPartition & domain ) override;

  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;

  /// the max number of fluid phases
  integer m_numPhases;

  /// the number of fluid components
  integer m_numComponents;

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// flag indicating whether total mass equation should be used
  integer m_useTotalMassEquation;

  /// list of relative permeability model names per target region
  array1d< string > m_relPermModelNames;

  /// maximum (absolute) change in a component fraction between two Newton iterations
  real64 m_maxCompFracChange;

  /// maximum (relative) change in pressure between two Newton iterations
  real64 m_maxRelativePresChange;

  /// maximum (absolute) change in pressure between two Newton iterations
  real64 m_maxAbsolutePresChange;

  /// maximum (relative) change in component density between two Newton iterations
  real64 m_maxRelativeCompDensChange;

  /// minimum value of the scaling factor obtained by enforcing maxCompFracChange
  real64 m_minScalingFactor;

  /// flag indicating whether local (cell-wise) chopping of negative compositions is allowed
  integer m_allowCompDensChopping;

  /// index of the target phase, used to impose the phase rate constraint
  localIndex m_targetPhaseIndex;

  /// name of the fluid constitutive model used as a reference for component/phase description
  string m_referenceFluidModelName;

};

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_
