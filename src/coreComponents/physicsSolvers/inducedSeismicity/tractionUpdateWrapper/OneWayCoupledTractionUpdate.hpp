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
 * @file OneWayCoupledTractionUpdate.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_ONEWAYCOUPLEDTRACTIONUPDATE_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_ONEWAYCOUPLEDTRACTIONUPDATE_HPP


#include "physicsSolvers/inducedSeismicity/tractionUpdateWrapper/FaultTractionUpdate.hpp"
#include "physicsSolvers/contact/ContactSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"


namespace geos
{

namespace inducedSeismicity
{

class OneWayCoupledTractionUpdate : public FaultTractionUpdate< ContactSolverBase, FlowSolverBase >
{
public:

  using Base = FaultTractionUpdate< ContactSolverBase, FlowSolverBase >;


  OneWayCoupledTractionUpdate( ContactSolverBase * contactSolver,
                               FlowSolverBase * flowSolver ):
    Base( contactSolver, flowSolver )
  {}

  virtual ~OneWayCoupledTractionUpdate() = default;

  virtual real64 updateFaultTraction( real64 const & time_n,
                                      real64 const & dt,
                                      const int cycleNumber,
                                      DomainPartition & domain ) const override final;

private:

  ContactSolverBase * getContactSolver() const { return std::get< 0 >( m_solvers ); }

  FlowSolverBase * getFlowSolver() const { return std::get< 1 >( m_solvers ); }

};

} // namespace inducedSeismicity

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_ONEWAYCOUPLEDTRACTIONUPDATE_HPP */
