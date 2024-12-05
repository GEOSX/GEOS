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
 * @file QuasiDynamicEQBase.cpp
 */

#include "QuasiDynamicEQBase.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mesh/DomainPartition.hpp"
#include "kernels/RateAndStateKernels.hpp"
#include "rateAndStateFields.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"

/// THIS is an alternative implementation to avoid the use of the TractionUpdateWrapper

namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

QuasiDynamicEQBase::QuasiDynamicEQBase( const string & name,
                                        Group * const parent ):
  PhysicsSolverBase( name, parent ),
  m_shearImpedance( 0.0 ),
  m_targetSlipIncrement( 1.0e-7 )
{
  this->registerWrapper( viewKeyStruct::shearImpedanceString(), &m_shearImpedance ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear impedance." );

  this->registerWrapper( viewKeyStruct::targetSlipIncrementString(), &m_targetSlipIncrement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0e-7 ).
    setDescription( "Target slip incrmeent for timestep size selction" );
}

QuasiDynamicEQBase::~QuasiDynamicEQBase()
{
  // TODO Auto-generated destructor stub
}

void QuasiDynamicEQBase::registerDataOnMesh( Group & meshBodies )
{
  PhysicsSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      SurfaceElementSubRegion & subRegion )
    {
      // Scalar functions on fault
      subRegion.registerField< rateAndState::stateVariable >( getName() );
      subRegion.registerField< rateAndState::stateVariable_n >( getName() );
      subRegion.registerField< rateAndState::slipRate >( getName() );
      subRegion.registerField< rateAndState::slipRate_n >( getName() );

      // Tangent (2-component) functions on fault
      string const labels2Comp[2] = {"tangent1", "tangent2" };
      subRegion.registerField< rateAndState::slipVelocity >( getName() ).
        setDimLabels( 1, labels2Comp ).reference().resizeDimension< 1 >( 2 );
      subRegion.registerField< rateAndState::shearTraction >( getName() ).
        setDimLabels( 1, labels2Comp ).reference().resizeDimension< 1 >( 2 );
    } );
  } );
}

void QuasiDynamicEQBase::applyInitialConditionsToFault( int const cycleNumber,
                                                        DomainPartition & domain ) const
{
  if( cycleNumber == 0 )
  {
    /// Apply initial conditions to the Fault
    FieldSpecificationManager & fieldSpecificationManager = FieldSpecificationManager::getInstance();

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )

    {
      fieldSpecificationManager.applyInitialConditions( mesh );
      mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                             [&]( localIndex const,
                                                                                  SurfaceElementSubRegion & subRegion )
      {
        arrayView1d< real64 const > const slipRate        = subRegion.getField< rateAndState::slipRate >();
        arrayView1d< real64 const > const stateVariable   = subRegion.getField< rateAndState::stateVariable >();
        arrayView1d< real64 > const stateVariable_n       = subRegion.getField< rateAndState::stateVariable_n >();
        arrayView1d< real64 > const slipRate_n            = subRegion.getField< rateAndState::slipRate_n >();

        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {
          slipRate_n [k]     = slipRate[k];
          stateVariable_n[k] = stateVariable[k];
        } );
      } );
    } );
  } 
}

void QuasiDynamicEQBase::solveRateAndStateEquations( real64 const time_n,
                                                     real64 const dt,
                                                     DomainPartition & domain ) const
{
  integer const maxNewtonIter = m_nonlinearSolverParameters.m_maxIterNewton;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion & subRegion )
    {
      // solve rate and state equations.
      rateAndStateKernels::createAndLaunch< parallelDevicePolicy<> >( subRegion, viewKeyStruct::frictionLawNameString(), m_shearImpedance, maxNewtonIter, time_n, dt );
      // save old state
      updateSlip( subRegion, dt );
    } );
  } );
}

void QuasiDynamicEQBase::updateSlip( ElementSubRegionBase & subRegion, real64 const dt ) const
{
  arrayView2d< real64 const > const slipVelocity    = subRegion.getField< rateAndState::slipVelocity >();
  arrayView2d< real64 > const deltaSlip             = subRegion.getField< contact::deltaSlip >();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    deltaSlip[k][0]     = slipVelocity[k][0] * dt;
    deltaSlip[k][1]     = slipVelocity[k][1] * dt;
  } );
}

real64 QuasiDynamicEQBase::setNextDt( real64 const & currentDt, DomainPartition & domain )
{
  GEOS_UNUSED_VAR( currentDt );

  real64 maxSlipRate = 0.0;
  // Spring-slider shear traction computation
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
