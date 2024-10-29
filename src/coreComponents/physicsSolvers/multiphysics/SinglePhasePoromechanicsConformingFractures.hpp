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
 * @file SinglePhasePoromechanicsConformingFractures.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_

#include "physicsSolvers/multiphysics/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/contact/SolidMechanicsLagrangeContact.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "dataRepository/Group.hpp"

namespace geos
{

template< typename FLOW_SOLVER = SinglePhaseBase >
class SinglePhasePoromechanicsConformingFractures : public SinglePhasePoromechanics< FLOW_SOLVER, SolidMechanicsLagrangeContact >
{
public:

  using Base = SinglePhasePoromechanics< FLOW_SOLVER, SolidMechanicsLagrangeContact >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "poromechanicsConformingFractures"; }

  /**
   * @brief main constructor for SinglePhasePoromechanicsConformingFractures objects
   * @param name the name of this instantiation of SinglePhasePoromechanicsConformingFractures in the repository
   * @param parent the parent group of this instantiation of SinglePhasePoromechanicsConformingFractures
   */
  SinglePhasePoromechanicsConformingFractures( const string & name,
                                               dataRepository::Group * const parent );

  /// Destructor for the class
  ~SinglePhasePoromechanicsConformingFractures() override {}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new SinglePhasePoromechanicsConformingFractures object through the object
   * catalog.
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< FLOW_SOLVER, SinglePhaseBase > )
    {
      return "SinglePhasePoromechanicsConformingFractures";
    }
    else
    {
      return FLOW_SOLVER::catalogName() + "PoromechanicsConformingFractures";
    }
  }

  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override final;

  virtual void updateState( DomainPartition & domain ) override final;

  virtual void setMGRStrategy() override
  {
    if( this->m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::mgr )
      GEOS_ERROR( GEOS_FMT( "{}: MGR strategy is not implemented for {}", this->getName(), this->getCatalogName()));
  }

  /**@}*/

private:

  struct viewKeyStruct : public Base::viewKeyStruct
  {};

  static const localIndex m_maxFaceNodes=11; // Maximum number of nodes on a contact face

  void assembleElementBasedContributions( real64 const time_n,
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs );

  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override final;

  void assembleForceResidualDerivativeWrtPressure( MeshLevel const & mesh,
                                                   arrayView1d< string const > const & regionNames,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs );

  void assembleFluidMassResidualDerivativeWrtDisplacement( MeshLevel const & mesh,
                                                           arrayView1d< string const > const & regionNames,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs );

  /**
   * @Brief add the nnz induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param rowLenghts the nnz in each row
   */
  void addTransmissibilityCouplingNNZ( DomainPartition const & domain,
                                       DofManager const & dofManager,
                                       arrayView1d< localIndex > const & rowLengths ) const;

  /**
   * @Brief add the sparsity pattern induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  void addTransmissibilityCouplingPattern( DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           SparsityPatternView< globalIndex > const & pattern ) const;

  /**
   * @brief Set up the Dflux_dApertureMatrix object
   *
   * @param domain
   * @param dofManager
   * @param localMatrix
   */
  void setUpDflux_dApertureMatrix( DomainPartition & domain,
                                   DofManager const & dofManager,
                                   CRSMatrix< real64, globalIndex > & localMatrix );

  /**
   * @brief
   *
   * @param domain
   */
  void updateHydraulicApertureAndFracturePermeability( DomainPartition & domain );

  std::unique_ptr< CRSMatrix< real64, localIndex > > & getRefDerivativeFluxResidual_dAperture()
  {
    return m_derivativeFluxResidual_dAperture;
  }

  CRSMatrixView< real64, localIndex const > getDerivativeFluxResidual_dNormalJump()
  {
    return m_derivativeFluxResidual_dAperture->toViewConstSizes();
  }

  CRSMatrixView< real64 const, localIndex const > getDerivativeFluxResidual_dNormalJump() const
  {
    return m_derivativeFluxResidual_dAperture->toViewConst();
  }

  std::unique_ptr< CRSMatrix< real64, localIndex > > m_derivativeFluxResidual_dAperture;

  string const m_pressureKey = SinglePhaseBase::viewKeyStruct::elemDofFieldString();
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_ */
