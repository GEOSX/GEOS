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



#include "OneWayCoupledTractionUpdate.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace inducedSeismicity
{


  real64 OneWayCoupledTractionUpdate::updateFaultTraction( real64 const & time_n,
                                                           real64 const & dt,
                                                           const int cycleNumber,
                                                           DomainPartition & domain ) const 
  {

    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->solverStep( time_n, dt, cycleNumber, domain );
    } );

    /// Now update the fault traction
    getContactSolver()->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                                     MeshLevel & mesh,
                                                                                     arrayView1d< string const > const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames, [&]( localIndex const, SurfaceElementSubRegion & subRegion )
      {
        arrayView2d< real64 > traction = subRegion.getField< fields::contact::traction >();
        arrayView1d< real64 > pressure = subRegion.getField< fields::flow::pressure >();
        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
        {
          // subtract pressure from the background normal stress
          traction( i, 0 ) -= pressure[i];
        } );
      } );
    } );
    return dt;
  }


} // namespace inducedSeismicity

} // namespace geos
