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
 * @file OneWayCoupledFractureFlowContactMechanics.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_ONEWAYCOUPLEDFRACTUREFLOWCONTACTMECHANICS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_ONEWAYCOUPLEDFRACTUREFLOWCONTACTMECHANICS_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/solidMechanics/contact/SolidMechanicsLagrangeContactBubbleStab.hpp"
#include "dataRepository/Group.hpp"

namespace geos
{

template< typename FLOW_SOLVER = SinglePhaseBase >
class OneWayCoupledFractureFlowContactMechanics : public CoupledSolver< FLOW_SOLVER, SolidMechanicsLagrangeContactBubbleStab >
{
public:

  using Base = CoupledSolver< FLOW_SOLVER, SolidMechanicsLagrangeContactBubbleStab >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  enum class SolverType : integer
  {
    Flow = 0,
    SolidMechanics = 1
  };

  /**
   * @brief main constructor for OneWayCoupledFractureFlowContactMechanics objects
   * @param name the name of this instantiation of OneWayCoupledFractureFlowContactMechanics in the repository
   * @param parent the parent group of this instantiation of OneWayCoupledFractureFlowContactMechanics
   */
  OneWayCoupledFractureFlowContactMechanics( const string & name,
                                             dataRepository::Group * const parent );

  /// Destructor for the class
  ~OneWayCoupledFractureFlowContactMechanics() override {}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new OneWayCoupledFractureFlowContactMechanics object through the object
   * catalog.
   */
  static string catalogName()
  {
    return "OneWayCoupledFractureFlowContactMechanics";
  }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual real64 sequentiallyCoupledSolverStep( real64 const & time_n,
                                                real64 const & dt,
                                                int const cycleNumber,
                                                DomainPartition & domain ) override final;

  virtual void postInputInitialization() override final;

  /**@}*/

private:

  struct viewKeyStruct : public Base::viewKeyStruct
  {};

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_ONEWAYCOUPLEDFRACTUREFLOWCONTACTMECHANICS_HPP_ */
