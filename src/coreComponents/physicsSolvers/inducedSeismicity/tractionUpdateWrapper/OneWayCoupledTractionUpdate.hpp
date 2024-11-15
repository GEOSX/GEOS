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

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_FAULTTRACTIONUPDATE_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_FAULTTRACTIONUPDATE_HPP


#include "physicsSolvers/inducedSeismicity/tractionUpdateWrapper/FaultTractionUpdate.hpp"


namespace geos
{

namespace inducedSeismicity
{  

template< typename CONTACT_SOLVER_TYPE , typename FLOW_SOLVER_TYPE >
class OneWayCoupledTractionUpdate : FaultTractionUpdate< CONTACT_SOLVER_TYPE, FLOW_SOLVER_TYPE>
{
public:

using Base = FaultTractionUpdate< CONTACT_SOLVER_TYPE, FLOW_SOLVER_TYPE>;


OneWayCoupledTractionUpdate( CONTACT_SOLVER_TYPE * contactSolver, 
                             FLOW_SOLVER_TYPE * flowSolver ) :
 Base( contactSolver, flowSolver )
{}

virtual ~OneWayCoupledTractionUpdate() = default;

virtual void registerMissingDataOnMesh() const = 0;

virtual void updateFaultTraction( real64 const & time_n,
                                  real64 const & dt,
                                  const int cycleNumber,
                                  DomainPartition & domain ) const override final
{
 getFlowSolver()->solverStep( time_n, dt, cycleNumber, domain );

 getContactSolver()->solverStep( time_n, dt, cycleNumber, domain );

 /// Now update the fault traction

 getContactSolver()-> forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
  {
    

  } );
}

private:

CONTACT_SOLVER_TYPE * getContactSolver() const { return std::get<0>( m_solvers ); }

FLOW_SOLVER_TYPE * getFlowSolver() const { return std::get<1>( m_solvers ); }

};

}

};

}

}

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_FAULTTRACTIONUPDATE_HPP */