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
 * @file QuasiDynamicEarthQuake.cpp
 */

#include "QuasiDynamicEarthQuake.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mesh/DomainPartition.hpp"
#include "rateAndStateFields.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"

#include "ExplicitQDRateAndState.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;
using namespace rateAndStateKernels;

template< typename RSSOLVER_TYPE >
QuasiDynamicEarthQuake< RSSOLVER_TYPE >::QuasiDynamicEarthQuake( const string & name,
                                                                 Group * const parent ):
  RSSOLVER_TYPE( name, parent ),
  m_stressSolverName(),
  m_stressSolver( nullptr )
{
  this->registerWrapper( viewKeyStruct::stressSolverNameString(), &m_stressSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of solver for computing stress." );
}

template< typename RSSOLVER_TYPE >
void QuasiDynamicEarthQuake< RSSOLVER_TYPE >::postInputInitialization()
{

  // Initialize member stress solver as specified in XML input
  m_stressSolver = &this->getParent().template getGroup< PhysicsSolverBase >( m_stressSolverName );

  PhysicsSolverBase::postInputInitialization();
}

template< typename RSSOLVER_TYPE >
QuasiDynamicEarthQuake< RSSOLVER_TYPE >::~QuasiDynamicEarthQuake()
{
  // TODO Auto-generated destructor stub
}

template< typename RSSOLVER_TYPE >
real64 QuasiDynamicEarthQuake< RSSOLVER_TYPE >::updateStresses( real64 const & time_n,
                                                                real64 const & dt,
                                                                const int cycleNumber,
                                                                DomainPartition & domain ) const
{
  // 1. Solve the momentum balance
  real64 const dtAccepted = m_stressSolver->solverStep( time_n, dt, cycleNumber, domain );
  // 2. Add background stress and possible forcing. 
  return dtAccepted;
}

// void QuasiDynamicEarthQuake::updateSlip( ElementSubRegionBase & subRegion, real64 const dt ) const
// {
//   arrayView2d< real64 const > const slipVelocity    = subRegion.getField< rateAndState::slipVelocity >();
//   arrayView2d< real64 > const deltaSlip             = subRegion.getField< contact::deltaSlip >();

//   arrayView2d< real64 > const dispJump = subRegion.getField< contact::targetIncrementalJump >();

//   forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
//   {
//     deltaSlip[k][0] = slipVelocity[k][0] * dt;
//     deltaSlip[k][1] = slipVelocity[k][1] * dt;
//     // Update tangential components of the displacement jump
//     dispJump[k][1] = slipVelocity[k][0] * dt;
//     dispJump[k][2] = slipVelocity[k][1] * dt;
//   } );
// }
template class QuasiDynamicEarthQuake< ImplicitQDRateAndState >;
template class QuasiDynamicEarthQuake< ExplicitQDRateAndState >;

namespace
{
typedef QuasiDynamicEarthQuake< ImplicitQDRateAndState > ImplicitQuasiDynamicEarthQuake;
typedef QuasiDynamicEarthQuake< ExplicitQDRateAndState > ExplicitQuasiDynamicEarthQuake;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, ImplicitQuasiDynamicEarthQuake, string const &, dataRepository::Group * const )
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, ExplicitQuasiDynamicEarthQuake, string const &, dataRepository::Group * const )
}

} // namespace geos