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

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_ONEWAYCOUPLEDTRACTIONUPDATE_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_ONEWAYCOUPLEDTRACTIONUPDATE_HPP


#include "physicsSolvers/inducedSeismicity/tractionUpdateWrapper/FaultTractionUpdate.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace inducedSeismicity
{  

template< >
class OneWayCoupledTractionUpdate : FaultTractionUpdate< ContactSolverBase, FlowSolverBase>
{
public:

using Base = FaultTractionUpdate< ContactSolverBase, FlowSolverBase>;


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

  forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
  {
    solver->solverStep( time_n, dt, cycleNumber, domain );
  } );

 /// Now update the fault traction
 getContactSolver()-> forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                      MeshLevel & mesh,
                                                      arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames, [&]( localIndex const, SurfaceElementSubRegion const & subRegion )
    {
      arrayView2d< real64 > traction = subregion.getField< fields::contact::traction >();
      arrayView1d< real64 > pressure = subregion.getField< fields::flow::pressure >();
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
      {
        // subtract pressure from the background normal stress
        traction( i, 0 ) = backgroundNormalStress - pressure( i );
      } );
    } );
}

private:

ContactSolverBase * getContactSolver() const { return std::get<0>( m_solvers ); }

FlowSolverBase * getFlowSolver() const { return std::get<1>( m_solvers ); }

};

} // namespace inducedSeismicity

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_ONEWAYCOUPLEDTRACTIONUPDATE_HPP */