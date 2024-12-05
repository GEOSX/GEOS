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
 * @file QuasiDynamicEQ.cpp
 */

#include "QuasiDynamicEQ.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mesh/DomainPartition.hpp"
#include "kernels/RateAndStateKernels.hpp"
#include "rateAndStateFields.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/inducedSeismicity/tractionUpdateWrapper/TractionUpdateFactory.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

QuasiDynamicEQ::QuasiDynamicEQ( const string & name,
                                Group * const parent ):
  PhysicsSolverBase( name, parent ),
  m_shearImpedance( 0.0 ),
  m_targetSlipIncrement( 1.0e-7 )
{
  this->registerWrapper( viewKeyStruct::shearImpedanceString(), &m_shearImpedance ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear impedance." );

  this->registerWrapper( viewKeyStruct::tractionUpdateTypeString(), &m_tractionUpdateType ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "." );

  this->registerWrapper( viewKeyStruct::contactSolverNameString(), &m_contactSolverName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "." );

  this->registerWrapper( viewKeyStruct::flowSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "." );

  this->registerWrapper( viewKeyStruct::targetSlipIncrementString(), &m_targetSlipIncrement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0e-7 ).
    setDescription( "Target slip incrmeent for timestep size selction" );
}

void QuasiDynamicEQ::postInputInitialization()
{

  // Initialize member stress solver as specified in XML input
  m_tractionUpdate = inducedSeismicity::tractionUpdateFactory( m_tractionUpdateType,
                                                               m_contactSolverName,
                                                               m_flowSolverName,
                                                               this );

  PhysicsSolverBase::postInputInitialization();
}

QuasiDynamicEQ::~QuasiDynamicEQ()
{}

void QuasiDynamicEQ::registerDataOnMesh( Group & meshBodies )
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

      m_tractionUpdate->registerMissingDataOnMesh( subRegion, this->getName() );
    } );
  } );
}

real64 QuasiDynamicEQ::solverStep( real64 const & time_n,
                                   real64 const & dt,
                                   int const cycleNumber,
                                   DomainPartition & domain )
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

  /// 1. Compute shear and normal tractions
  GEOS_LOG_LEVEL_RANK_0( 1, "Stress solver" );

  real64 const dtStress = m_tractionUpdate->updateFaultTraction( time_n, dt, cycleNumber, domain );

  /// 2. Solve for slip rate and state variable and, compute slip
  GEOS_LOG_LEVEL_RANK_0( 1, "Rate and State solver" );

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
      rateAndStateKernels::createAndLaunch< parallelDevicePolicy<> >( subRegion, viewKeyStruct::frictionLawNameString(), m_shearImpedance, maxNewtonIter, time_n, dtStress );
      // save old state
      updateSlip( subRegion, dtStress );
    } );
  } );

  // return time step size achieved by stress solver
  return dtStress;
}


void QuasiDynamicEQ::updateSlip( ElementSubRegionBase & subRegion, real64 const dt ) const
{
  arrayView2d< real64 const > const slipVelocity    = subRegion.getField< rateAndState::slipVelocity >();
  arrayView2d< real64 > const deltaSlip             = subRegion.getField< contact::deltaSlip >();

  arrayView2d< real64 > const dispJump = subRegion.getField< contact::targetIncrementalJump >();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    deltaSlip[k][0]     = slipVelocity[k][0] * dt;
    deltaSlip[k][1]     = slipVelocity[k][1] * dt;
    // Update tangential components of the displacement jump
    dispJump[k][1]      = slipVelocity[k][0] * dt;
    dispJump[k][2]      = slipVelocity[k][1] * dt;
  } );
}

real64 QuasiDynamicEQ::setNextDt( real64 const & currentDt, DomainPartition & domain )
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

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, QuasiDynamicEQ, string const &, dataRepository::Group * const )

} // namespace geos
