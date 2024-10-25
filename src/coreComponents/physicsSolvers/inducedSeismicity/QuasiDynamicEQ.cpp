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

namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

QuasiDynamicEQ::QuasiDynamicEQ( const string & name,
                                Group * const parent ):
  SolverBase( name, parent ),
  m_stressSolver( nullptr ),
  m_stressSolverName( "SpringSlider" ),
  m_shearImpedance( 0.0 ),
  m_targetSlipIncrement( 1.0e-7 )
{
  this->registerWrapper( viewKeyStruct::shearImpedanceString(), &m_shearImpedance ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear impedance." );

  this->registerWrapper( viewKeyStruct::stressSolverNameString(), &m_stressSolverName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of solver for computing stress. If empty, the spring-slider model is run." );

  this->registerWrapper( viewKeyStruct::targetSlipIncrementString(), &m_targetSlipIncrement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0e-7 ).
    setDescription( "Target slip incrmeent for timestep size selction" );
}

void QuasiDynamicEQ::postInputInitialization()
{

  // Initialize member stress solver as specified in XML input
  if( !m_stressSolverName.empty() )
  {
    m_stressSolver = &this->getParent().getGroup< SolverBase >( m_stressSolverName );
  }

  SolverBase::postInputInitialization();
}

QuasiDynamicEQ::~QuasiDynamicEQ()
{
  // TODO Auto-generated destructor stub
}

void QuasiDynamicEQ::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

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

      // Tangent (2-component) functions on fault
      string const labels2Comp[2] = {"tangent1", "tangent2" };
      subRegion.registerField< rateAndState::slipVelocity >( getName() ).
        setDimLabels( 1, labels2Comp ).reference().resizeDimension< 1 >( 2 );
      subRegion.registerField< rateAndState::deltaSlip >( getName() ).
        setDimLabels( 1, labels2Comp ).reference().resizeDimension< 1 >( 2 );

      if( !subRegion.hasWrapper( contact::dispJump::key() ))
      {
        // 3-component functions on fault
        string const labels3Comp[3] = { "normal", "tangent1", "tangent2" };
        subRegion.registerField< contact::dispJump >( getName() ).
          setDimLabels( 1, labels3Comp ).
          reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< contact::traction >( getName() ).
          setDimLabels( 1, labels3Comp ).
          reference().resizeDimension< 1 >( 3 );

        subRegion.registerWrapper< string >( viewKeyStruct::frictionLawNameString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRestartFlags( RestartFlags::NO_WRITE ).
          setSizedFromParent( 0 );

        string & frictionLawName = subRegion.getReference< string >( viewKeyStruct::frictionLawNameString() );
        frictionLawName = SolverBase::getConstitutiveName< FrictionBase >( subRegion );
        GEOS_ERROR_IF( frictionLawName.empty(), GEOS_FMT( "{}: FrictionBase model not found on subregion {}",
                                                          getDataContext(), subRegion.getDataContext() ) );
      }
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
                                                                 arrayView1d< string const > const & )

    {
      fieldSpecificationManager.applyInitialConditions( mesh );
    } );
  }

  /// 1. Compute shear and normal tractions
  GEOS_LOG_LEVEL_RANK_0( 1, "Stress solver" );

  real64 const dtStress = updateStresses( time_n, dt, cycleNumber, domain );

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
      saveOldStateAndUpdateSlip( subRegion, dtStress );
    } );
  } );

  // return time step size achieved by stress solver
  return dtStress;
}

real64 QuasiDynamicEQ::updateStresses( real64 const & time_n,
                                       real64 const & dt,
                                       const int cycleNumber,
                                       DomainPartition & domain ) const
{
  // Call member variable stress solver to update the stress state
  if( m_stressSolver )
  {
    // 1. Solve the momentum balance
    real64 const dtStress =  m_stressSolver->solverStep( time_n, dt, cycleNumber, domain );

    return dtStress;
  }
  else
  {
    // Spring-slider shear traction computation
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )

    {
      mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                             [&]( localIndex const,
                                                                                  SurfaceElementSubRegion & subRegion )
      {

        arrayView2d< real64 const > const deltaSlip = subRegion.getField< rateAndState::deltaSlip >();
        arrayView2d< real64 > const traction        = subRegion.getField< fields::contact::traction >();

        string const & fricitonLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
        RateAndStateFriction const & frictionLaw = getConstitutiveModel< RateAndStateFriction >( subRegion, fricitonLawName );

        RateAndStateFriction::KernelWrapper frictionKernelWrapper = frictionLaw.createKernelUpdates();

        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {
          SpringSliderParameters springSliderParameters = SpringSliderParameters( traction[k][0],
                                                                                  frictionKernelWrapper.getACoefficient( k ),
                                                                                  frictionKernelWrapper.getBCoefficient( k ),
                                                                                  frictionKernelWrapper.getDcCoefficient( k ) );


          traction[k][1] = traction[k][1] + springSliderParameters.tauRate * dt
                           - springSliderParameters.springStiffness * deltaSlip[k][0];
          traction[k][2] = traction[k][2] + springSliderParameters.tauRate * dt
                           - springSliderParameters.springStiffness * deltaSlip[k][1];
        } );
      } );
    } );
    return dt;
  }
}

void QuasiDynamicEQ::saveOldStateAndUpdateSlip( ElementSubRegionBase & subRegion, real64 const dt ) const
{
  arrayView1d< real64 > const stateVariable   = subRegion.getField< rateAndState::stateVariable >();
  arrayView1d< real64 > const stateVariable_n = subRegion.getField< rateAndState::stateVariable_n >();
  arrayView2d< real64 > const slipVelocity    = subRegion.getField< rateAndState::slipVelocity >();
  arrayView2d< real64 > const deltaSlip       = subRegion.getField< rateAndState::deltaSlip >();

  arrayView2d< real64 > const dispJump = subRegion.getField< contact::dispJump >();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    stateVariable_n[k]  = stateVariable[k];
    deltaSlip[k][0]     = slipVelocity[k][0] * dt;
    deltaSlip[k][1]     = slipVelocity[k][1] * dt;
    // Update tangential components of the displacement jump
    dispJump[k][1]      = dispJump[k][1] + slipVelocity[k][0] * dt;
    dispJump[k][2]      = dispJump[k][2] + slipVelocity[k][1] * dt;
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

REGISTER_CATALOG_ENTRY( SolverBase, QuasiDynamicEQ, string const &, dataRepository::Group * const )

} // namespace geos
