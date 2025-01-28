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
 * @file SinglePhasePoromechanicsConformingFracturesALM.cpp
 */

#include "SinglePhasePoromechanicsConformingFracturesALM.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

template< typename FLOW_SOLVER >
SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::SinglePhasePoromechanicsConformingFracturesALM( const string & name,
                                                                                                               Group * const parent )
  : Base( name, parent )
{}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::setupCoupling( DomainPartition const & domain,
                                                                                   DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( domain, dofManager );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::setupSystem( DomainPartition & domain,
                                                                                 DofManager & dofManager,
                                                                                 CRSMatrix< real64, globalIndex > & localMatrix,
                                                                                 ParallelVector & rhs,
                                                                                 ParallelVector & solution,
                                                                                 bool const setSparsity )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( domain, dofManager, localMatrix, rhs, solution, setSparsity );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::assembleSystem( real64 const time_n,
                                                                                 real64 const dt,
                                                                                 DomainPartition & domain,
                                                                                 DofManager const & dofManager,
                                                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                 arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n, dt, domain, dofManager, localMatrix, localRhs );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFracturesALM< FLOW_SOLVER >::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( domain );
}

template class SinglePhasePoromechanicsConformingFracturesALM<>;
template class SinglePhasePoromechanicsConformingFracturesALM< SinglePhaseReservoirAndWells<> >;

namespace
{
typedef SinglePhasePoromechanicsConformingFracturesALM< SinglePhaseReservoirAndWells<> > SinglePhaseReservoirPoromechanicsConformingFracturesALM;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseReservoirPoromechanicsConformingFracturesALM, string const &, Group * const )
typedef SinglePhasePoromechanicsConformingFracturesALM<> SinglePhasePoromechanicsConformingFracturesALM;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhasePoromechanicsConformingFracturesALM, string const &, Group * const )
}

} /* namespace geos */
