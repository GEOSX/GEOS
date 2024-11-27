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

/// THIS is an alternative implementation to avoid the use of the TractionUpdateWrapper

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEQBASE_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEQBASE_HPP

#include "physicsSolvers/PhysicsSolverBase.hpp"

namespace geos
{

template< typename RESERVOIR_SOLVER >
class CoupledReservoirEQ : public CoupledSolver< RESERVOIR_SOLVER, QuasiDynamicEQBase >
{
public:

  using Base = CoupledReservoirEQ< RESERVOIR_SOLVER, QuasiDynamicEQBase >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  enum class SolverType : integer
  {
    Reservoir = 0,
    Earthquake = 1
  };

  static string coupledSolverAttributePrefix() { return "Coupled"; }

  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  CoupledReservoirEQ() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  CoupledReservoirEQ( const string & name,
                      Group * const parent );

  /// Destructor
  virtual ~CoupledReservoirEQ() override;

  static string catalogName() { return "CoupledReservoirEQ"; }

  /**
   * @return Get the final class Catalog name
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override final

  {
    /// 1. Compute shear and normal tractions
    GEOS_LOG_LEVEL_RANK_0( 1, "Stress solver" );

    real64 const dtStress = reservoirSolver()->solverStep( time_n, dt, cycleNumber, domain );

    /// 2. Solve for slip rate and state variable and, compute slip
    GEOS_LOG_LEVEL_RANK_0( 1, "Earthquake solver" );
    earthquakeSolver()->solverStep( time_n, dt, cycleNumber, domain );
  }

  virtual real64 setNextDt( real64 const & currentDt,
                            DomainPartition & domain ) override final

  {
    return earthquakeSolver()->setNextDt( currentDt, domain );
  }

  /**
   * @brief accessor for the pointer to the reservoir solver
   * @return a pointer to the reservoir solver
   */
  RESERVOIR_SOLVER * reservoirSolver() const
  {
    return std::get< toUnderlying( SolverType::Reservoir ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the earthquake solver
   * @return a pointer to the earthquake solver
   */
  QuasiDynamicEQBase * earthquakeSolver() const
  {
    return std::get< toUnderlying( SolverType::Earthquake ) >( m_solvers );
  }

private:

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEQBASE_HPP */
