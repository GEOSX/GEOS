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
 * @file FlowProppantTransportSolver.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_FLOWPROPPANTTRANSPORTSOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_FLOWPROPPANTTRANSPORTSOLVER_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"

namespace geos
{

class ProppantTransport;
class FlowSolverBase;

class FlowProppantTransportSolver : public CoupledSolver< ProppantTransport,
                                                          FlowSolverBase >
{
public:

  using Base = CoupledSolver< ProppantTransport, FlowSolverBase >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  enum class SolverType : integer
  {
    ProppantTransport = 0,
    Flow = 1
  };

  /**
   * @brief main constructor for FlowProppantTransportSolver Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of FlowProppantTransportSolver
   */
  FlowProppantTransportSolver( const string & name,
                               Group * const parent );

  /// Destructor for the class
  ~FlowProppantTransportSolver() override {};

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new FlowProppantTransportSolver object through the object catalog.
   */
  static string catalogName() { return "FlowProppantTransport"; }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /**
   * @brief accessor for the pointer to the proppant transport solver
   * @return a pointer to the proppant transport solver
   */
  ProppantTransport * proppantTransportSolver() const
  {
    return std::get< toUnderlying( SolverType::ProppantTransport ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the flow solver
   * @return a pointer to the flow solver
   */
  FlowSolverBase * flowSolver() const
  {
    return std::get< toUnderlying( SolverType::Flow ) >( m_solvers );
  }

  /**@}*/

private:

  real64 sequentiallyCoupledSolverStep( real64 const & time_n,
                                        real64 const & dt,
                                        int const cycleNumber,
                                        DomainPartition & domain ) override final;

  /**
   * @brief Utility function to perform the pre-step Update
   * @param[in] time_n the time at the previous converged time step
   * @param[in] dt the time step size
   * @param[in] domain the domain partition
   */
  void preStepUpdate( real64 const & time_n,
                      real64 const & dt,
                      DomainPartition & domain );

  /**
   * @brief Utility function to perform the post-step Update
   * @param[in] time_n the time at the previous converged time step
   * @param[in] dt the time step size
   * @param[in] domain the domain partition
   */
  void postStepUpdate( real64 const & time_n,
                       real64 const & dt,
                       DomainPartition & domain );

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_FLOWPROPPANTTRANSPORTSOLVER_HPP_ */
