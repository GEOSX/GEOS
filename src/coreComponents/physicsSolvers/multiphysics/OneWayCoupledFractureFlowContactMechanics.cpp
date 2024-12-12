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
 * @file OneWayCoupledFractureFlowContactMechanics.cpp
 */

#include "OneWayCoupledFractureFlowContactMechanics.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "mesh/DomainPartition.hpp"


namespace geos
{
using namespace dataRepository;

template< typename FLOW_SOLVER >
OneWayCoupledFractureFlowContactMechanics< FLOW_SOLVER >::OneWayCoupledFractureFlowContactMechanics( const string & name,
                                                                                                     Group * const parent )
: Base( name, parent )
{}

template< typename FLOW_SOLVER >
void OneWayCoupledFractureFlowContactMechanics< FLOW_SOLVER >::postInputInitialization() 
{
  bool const isSequential = getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::Sequential;
  GEOS_THROW_IF( !isSequential ,
                   GEOS_FMT( "Only a sequential coupling is allowed for this solver." ),
                   InputError );

  Base::postInputInitialization();                 
}


template< typename FLOW_SOLVER >
real64 OneWayCoupledFractureFlowContactMechanics< FLOW_SOLVER >::sequentiallyCoupledSolverStep( real64 const & time_n,
                                                                                                real64 const & dt,
                                                                                                int const cycleNumber,
                                                                                                DomainPartition & domain )
{
  forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
  {
    solver->solverStep( time_n, dt, cycleNumber, domain );
  } );
  
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                     MeshLevel & mesh,
                                                                     arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion & subRegion )
    {
      auto traction = subRegion.getField< fields::contact::traction >();
      auto pressure = subRegion.getField< fields::flow::pressure >();
      
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        traction( k, 0 ) = traction( k, 0 ) + pressure[k];
      } );
    } );
  } );

  return dt;
}

template class OneWayCoupledFractureFlowContactMechanics< SinglePhaseBase >;

namespace
{
typedef OneWayCoupledFractureFlowContactMechanics< SinglePhaseBase > OneWayCoupledFractureFlowContactMechanicsSinglePhase;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, OneWayCoupledFractureFlowContactMechanicsSinglePhase, string const &, dataRepository::Group * const )
}

} /* namespace geos */
