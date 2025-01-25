/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhasePoromechanicsConformingFracturesALM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_

#include "physicsSolvers/multiphysics/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/contact/SolidMechanicsAugmentedLagrangianContact.hpp"

namespace geos
{

template< typename FLOW_SOLVER = SinglePhaseBase >
class SinglePhasePoromechanicsConformingFracturesALM : public SinglePhasePoromechanics< FLOW_SOLVER, SolidMechanicsAugmentedLagrangianContact >
{
public:

  using Base = SinglePhasePoromechanics< FLOW_SOLVER, SolidMechanicsAugmentedLagrangianContact >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "poromechanicsConformingFracturesALM"; }

  /**
   * @brief main constructor for SinglePhasePoromechanicsConformingFracturesALM objects
   * @param name the name of this instantiation of SinglePhasePoromechanicsConformingFracturesALM in the repository
   * @param parent the parent group of this instantiation of SinglePhasePoromechanicsConformingFracturesALM
   */
  SinglePhasePoromechanicsConformingFracturesALM( const string & name,
                                               dataRepository::Group * const parent );

  /// Destructor for the class
  ~SinglePhasePoromechanicsConformingFracturesALM() override {}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new SinglePhasePoromechanicsConformingFracturesALM object through the object
   * catalog.
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< FLOW_SOLVER, SinglePhaseBase > )
    {
      return "SinglePhasePoromechanicsConformingFracturesALM";
    }
    else
    {
      return FLOW_SOLVER::catalogName() + "PoromechanicsConformingFracturesALM";
    }
  }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override final;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override final;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override final;

  virtual void updateState( DomainPartition & domain ) override final;

  virtual void setMGRStrategy() override final
  {
    if( this->m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::mgr )
      GEOS_ERROR( GEOS_FMT( "{}: MGR strategy is not implemented for {}", this->getName(), this->getCatalogName()));
  }

  /**@}*/

private:

  struct viewKeyStruct : public Base::viewKeyStruct
  {};

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_ */
