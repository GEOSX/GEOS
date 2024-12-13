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
 * @file ImplicitQDRateAndState.cpp
 */

#include "ImplicitQDRateAndState.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mesh/DomainPartition.hpp"
#include "kernels/ImplicitRateAndStateKernels.hpp"
#include "rateAndStateFields.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"


namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

ImplicitQDRateAndState::ImplicitQDRateAndState( const string & name,
                                                Group * const parent ):
  QDRateAndStateBase( name, parent ),
  m_targetSlipIncrement( 1.0e-7 )
{
  this->registerWrapper( viewKeyStruct::targetSlipIncrementString(), &m_targetSlipIncrement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0e-7 ).
    setDescription( "Target slip incrmeent for timestep size selction" );
}

ImplicitQDRateAndState::~ImplicitQDRateAndState()
{
  // TODO Auto-generated destructor stub
}

void ImplicitQDRateAndState::solveRateAndStateEquations( real64 const time_n,
                                                         real64 const dt,
                                                         DomainPartition & domain ) const
{
  integer const maxNewtonIter = m_nonlinearSolverParameters.m_maxIterNewton;
  real64 const newtonTol = m_nonlinearSolverParameters.m_newtonTol;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion & subRegion )
    {
      // solve rate and state equations.
      rateAndStateKernels::createAndLaunch< rateAndStateKernels::ImplicitFixedStressRateAndStateKernel, parallelDevicePolicy<> >( subRegion,
                                                                                                                                  viewKeyStruct::frictionLawNameString(),
                                                                                                                                  m_shearImpedance,
                                                                                                                                  maxNewtonIter, newtonTol,
                                                                                                                                  time_n,
                                                                                                                                  dt ); // save
                                                                                                                                        // old
                                                                                                                                        // state
      updateSlip( subRegion, dt );
    } );
  } );
}

real64 ImplicitQDRateAndState::solverStep( real64 const & time_n,
                                           real64 const & dt,
                                           int const cycleNumber,
                                           DomainPartition & domain )
{
  applyInitialConditionsToFault( cycleNumber, domain );
  GEOS_LOG_LEVEL_RANK_0( 1, "Stress solver" );
  updateStresses( time_n, dt, cycleNumber, domain );
  GEOS_LOG_LEVEL_RANK_0( 1, "Rate and state solver" );
  solveRateAndStateEquations( time_n, dt, domain );
  saveState( domain );
  return dt;
}

void ImplicitQDRateAndState::updateSlip( ElementSubRegionBase & subRegion, real64 const dt ) const
{
  arrayView2d< real64 const > const slipVelocity    = subRegion.getField< rateAndState::slipVelocity >();
  arrayView2d< real64 > const deltaSlip             = subRegion.getField< contact::deltaSlip >();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    deltaSlip[k][0] = slipVelocity[k][0] * dt;
    deltaSlip[k][1] = slipVelocity[k][1] * dt;
  } );
}

real64 ImplicitQDRateAndState::setNextDt( real64 const & currentDt, DomainPartition & domain )
{
  GEOS_UNUSED_VAR( currentDt );

  real64 maxSlipRate = 0.0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    real64 maxSlipRateOnThisRank  = 0.0;
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion const & subRegion )
    {
      arrayView1d< real64 const > const slipRate = subRegion.getField< rateAndState::slipRate >();

      RAJA::ReduceMax< parallelDeviceReduce, real64 > maximumSlipRateOnThisRegion( 0.0 );
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        maximumSlipRateOnThisRegion.max( slipRate[k] );
      } );
      if( maximumSlipRateOnThisRegion.get() > maxSlipRateOnThisRank )
        maxSlipRateOnThisRank = maximumSlipRateOnThisRegion.get();
    } );
    maxSlipRate = MpiWrapper::max( maxSlipRateOnThisRank );
  } );

  real64 const nextDt = m_targetSlipIncrement / maxSlipRate;

  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "The next dt will be {:.2e} s", nextDt ));

  return nextDt;
}

} // namespace geos
