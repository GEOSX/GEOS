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
 * @file QuasiDynamicEQRK32.cpp
 */

#include "QuasiDynamicEQRK32.hpp"

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
using namespace rateAndStateKernels;

QuasiDynamicEQRK32::QuasiDynamicEQRK32( const string & name,
                                        Group * const parent ):
  SolverBase( name, parent ),
  m_stressSolver( nullptr ),
  m_stressSolverName( "SpringSlider" ),
  m_shearImpedance( 0.0 ),
  m_timestepAbsTol( 1.0e-5 ),
  m_timestepRelTol( 1.0e-5 ),

  m_timestepAcceptSafety( 0.81 ),
  m_prevTimestepErrors{ 0.0, 0.0 } ,
  m_beta{ 1.0/18.0, 1.0/9.0, 1.0/18.0 } ,
  m_rkOrders{ 3, 2 },
  m_successfulStep( false )
{
  this->registerWrapper( viewKeyStruct::shearImpedanceString(), &m_shearImpedance ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear impedance." );

  this->registerWrapper( viewKeyStruct::stressSolverNameString(), &m_stressSolverName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of solver for computing stress. If empty, the spring-slider model is run." );
}

void QuasiDynamicEQRK32::postInputInitialization()
{

  // Initialize member stress solver as specified in XML input
  if( !m_stressSolverName.empty() )
  {
    m_stressSolver = &this->getParent().getGroup< SolverBase >( m_stressSolverName );
  }

  SolverBase::postInputInitialization();
}

QuasiDynamicEQRK32::~QuasiDynamicEQRK32()
{
  // TODO Auto-generated destructor stub
}


// TODO The vectors listed below are only temporary memory buffers used within an adaptive
// time step in the solverStep function. They don't need to persist between time steps,
// and could be allocated in the beginning of each call to solverStep, instead of 
// being managed by the region manager. This would of course result in more heap allocations.
// Not sure what is preferable. 
// 
// slipVelocity_n
// traction_n
// deltaSlip_n
// stateVariable_n
void QuasiDynamicEQRK32::registerDataOnMesh( Group & meshBodies )
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
      subRegion.registerField< rateAndState::slipVelocity_n >( getName() ).
        setDimLabels( 1, labels2Comp ).reference().resizeDimension< 1 >( 2 );  
      subRegion.registerField< rateAndState::deltaSlip >( getName() ).
        setDimLabels( 1, labels2Comp ).reference().resizeDimension< 1 >( 2 );
      subRegion.registerField< rateAndState::deltaSlip_n >( getName() ).
        setDimLabels( 1, labels2Comp ).reference().resizeDimension< 1 >( 2 );  

      // Runge-Kutta stage rates
      subRegion.registerField< rateAndState::stateVariableRKStageRate >( getName() ).reference().resizeDimension< 1 >( 2 ); 
      subRegion.registerField< rateAndState::deltaSlipRKStageRate >( getName() ).reference().resizeDimension< 1, 2>( 2, 2 );
      // Error
      string const labelsError[3] = { "deltaSlip1", "deltaSlip2", "stateVariable"};
      subRegion.registerField< rateAndState::error >( getName() ).
        setDimLabels( 1, labelsError ).reference().resizeDimension< 1 >( 3 ); 
      

      if( !subRegion.hasWrapper( contact::dispJump::key() ))
      {
        // 3-component functions on fault
        string const labels3Comp[3] = { "normal", "tangent1", "tangent2" };
        subRegion.registerField< contact::dispJump >( getName() ).
          setDimLabels( 1, labels3Comp ).
          reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< contact::dispJump_n >( getName() ).
          setDimLabels( 1, labels3Comp ).
          reference().resizeDimension< 1 >( 3 );  
        subRegion.registerField< contact::traction >( getName() ).
          setDimLabels( 1, labels3Comp ).
          reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< contact::traction_n >( getName() ).
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

real64 QuasiDynamicEQRK32::solverStep( real64 const & time_n,
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
    saveState(domain);
  }

  real64 dtAdaptive = dt;

  GEOS_LOG_LEVEL_RANK_0( 1, "Begin adaptive time step" );
  while (true) // Adaptive time step loop. Performs a Runge-Kutta time stepping with error control on state and slip
  {
    
    // First Runge-Kutta stage
    
    // Evolve deltaSlip, dispJump and stateVariable to t + 0.5*dtAdaptive
    // Store the first stage rates for deltaSlip and stateVariable
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )

    {
      mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                            [&]( localIndex const,
                                                                                  SurfaceElementSubRegion & subRegion )
      {
        
        string const & fricitonLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
        RateAndStateFriction const & frictionLaw = getConstitutiveModel< RateAndStateFriction >( subRegion, fricitonLawName );
        RateAndStateFriction::KernelWrapper frictionKernelWrapper = frictionLaw.createKernelUpdates();

        // Fields read for this stage
        arrayView2d< real64 const > const slipVelocity_n  = subRegion.getField< rateAndState::slipVelocity >();
        arrayView1d< real64 const > const stateVariable_n = subRegion.getField< rateAndState::stateVariable_n >();
        arrayView2d< real64 const > const deltaSlip_n     = subRegion.getField< rateAndState::deltaSlip_n >();
        arrayView2d< real64 const > const dispJump_n      = subRegion.getField< contact::dispJump_n >();
        // Fields written for this stage
        arrayView1d< real64 > const stateVariable         = subRegion.getField< rateAndState::stateVariable >();
        arrayView1d< real64 > const slipRate              = subRegion.getField< rateAndState::slipRate >();
        arrayView2d< real64 > const deltaSlip             = subRegion.getField< rateAndState::deltaSlip >();
        arrayView2d< real64 > const dispJump              = subRegion.getField< contact::dispJump >();
        // Stage rates
        arrayView2d< real64 > const stateVariableRKStageRate  = subRegion.getField< rateAndState::stateVariableRKStageRate >();
        arrayView3d< real64 > const deltaSlipRKStageRate      = subRegion.getField< rateAndState::deltaSlipRKStageRate >();
        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {   
            slipRate[k] = LvArray::math::sqrt(slipVelocity_n[k][0] * slipVelocity_n[k][0] + slipVelocity_n[k][1]*slipVelocity_n[k][1]);
            // Runge-Kutta stage rates for stage 1
            stateVariableRKStageRate[k][0] = frictionKernelWrapper.stateEvolution(k, slipRate[k], stateVariable_n[k]);
            LvArray::tensorOps::copy< 2 >( deltaSlipRKStageRate[k][0], slipVelocity_n[k] );
            // Runge-Kutta stage values for stage 1
            stateVariable[k] = stateVariable_n[k] + 0.5*dtAdaptive*stateVariableRKStageRate[0][k];
            deltaSlip[k][0]  = deltaSlip_n[k][0]  + 0.5*dtAdaptive*deltaSlipRKStageRate[k][0][0];
            deltaSlip[k][1]  = deltaSlip_n[k][1]  + 0.5*dtAdaptive*deltaSlipRKStageRate[k][0][1];
            // Set tangential components of the displacement jump for stress solver and state variable for non-linear solve
            dispJump[k][1] = dispJump_n[k][1] + deltaSlip[k][0];
            dispJump[k][2] = dispJump_n[k][2] + deltaSlip[k][1];
        } );
      } );
    } );

    // Update stresses using the updated value of deltaSlip/dispJump at t + 0.5*dtAdapt
    real64 dtStress = updateStresses( time_n, 0.5*dtAdaptive, cycleNumber, domain );
    GEOS_UNUSED_VAR( dtStress );
    
    // Update slipRate/slipVelocity using tractions, and stateVariable at t + 0.5*dtAdapt
    updateSlipVelocity( time_n, 0.5*dt, domain );


    // Second Runge-Kutta stage
    
    // Evolve deltaSlip, dispJump and stateVariable to t + dtAdaptive
    // Store the second stage rates for deltaSlip and stateVariable
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )

    {
      mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                            [&]( localIndex const,
                                                                                  SurfaceElementSubRegion & subRegion )
      {
        
        string const & fricitonLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
        RateAndStateFriction const & frictionLaw = getConstitutiveModel< RateAndStateFriction >( subRegion, fricitonLawName );
        RateAndStateFriction::KernelWrapper frictionKernelWrapper = frictionLaw.createKernelUpdates();

        // Fields read for this stage
        arrayView1d< real64 const > const stateVariable_n       = subRegion.getField< rateAndState::stateVariable_n >();
        arrayView2d< real64 const > const deltaSlip_n           = subRegion.getField< rateAndState::deltaSlip_n >();
        arrayView2d< real64 const > const dispJump_n            = subRegion.getField< contact::dispJump_n >();
        arrayView1d< real64 const > const slipRate              = subRegion.getField< rateAndState::slipRate >();
        arrayView2d< real64 const > const slipVelocity          = subRegion.getField< rateAndState::slipVelocity >();
        
        // Fields written for this stage
        arrayView1d< real64 > const stateVariable         = subRegion.getField< rateAndState::stateVariable >();
        arrayView2d< real64 > const deltaSlip             = subRegion.getField< rateAndState::deltaSlip >();
        arrayView2d< real64 > const dispJump              = subRegion.getField< contact::dispJump >();

        // Stage rates
        arrayView2d< real64 > const stateVariableRKStageRate  = subRegion.getField< rateAndState::stateVariableRKStageRate >();
        arrayView3d< real64 > const deltaSlipRKStageRate      = subRegion.getField< rateAndState::deltaSlipRKStageRate >();

        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {   
            // Runge-Kutta stage rates for stage 2
            stateVariableRKStageRate[k][1] = frictionKernelWrapper.stateEvolution(k, slipRate[k], stateVariable[k]);
            LvArray::tensorOps::copy< 2 >( deltaSlipRKStageRate[k][1], slipVelocity[k] );
            
            // Runge-Kutta stage values for stage 2
            stateVariable[k] = stateVariable_n[k] + dtAdaptive*( -stateVariableRKStageRate[k][0] + 2*stateVariableRKStageRate[k][1] );
            deltaSlip[k][0]  = deltaSlip_n[k][0]  + dtAdaptive*( -deltaSlipRKStageRate[k][0][0]  + 2*deltaSlipRKStageRate[k][1][0] );
            deltaSlip[k][1]  = deltaSlip_n[k][1]  + dtAdaptive*( -deltaSlipRKStageRate[k][0][1]  + 2*deltaSlipRKStageRate[k][1][1] );
            // Set tangential components of the displacement jump for stress solver and state variable for non-linear solve
            dispJump[k][1] = dispJump_n[k][1] + deltaSlip[k][0];
            dispJump[k][2] = dispJump_n[k][2] + deltaSlip[k][1];
        } );
      } );
    } );

    // Update stresses using the updated value of deltaSlip/dispJump at dtAdapt
    dtStress = updateStresses( time_n, dtAdaptive, cycleNumber, domain );
    
    // Update slipRate/slipVelocity using tractions, and stateVariable at t + 0.5*dtAdapt
    updateSlipVelocity( time_n, dtAdaptive, domain );

    // Third Runge-Kutta stage
    
    // Perform third and second order update of deltaSlip and stateVariable, and compute
    // error
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )

    {
      mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                            [&]( localIndex const,
                                                                                  SurfaceElementSubRegion & subRegion )
      {
        
        string const & fricitonLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
        RateAndStateFriction const & frictionLaw = getConstitutiveModel< RateAndStateFriction >( subRegion, fricitonLawName );
        RateAndStateFriction::KernelWrapper frictionKernelWrapper = frictionLaw.createKernelUpdates();
        
        // Fields read for this stage
        arrayView1d< real64 const > const stateVariable_n = subRegion.getField< rateAndState::stateVariable_n >();
        arrayView2d< real64 const > const deltaSlip_n     = subRegion.getField< rateAndState::deltaSlip_n >();
        arrayView1d< real64 const > const slipRate        = subRegion.getField< rateAndState::slipRate >();
        arrayView2d< real64 const > const slipVelocity    = subRegion.getField< rateAndState::slipVelocity >();
        arrayView2d< real64 const > const dispJump_n      = subRegion.getField< contact::dispJump_n >();
        
        // Fields written for this stage
        arrayView1d< real64 > const stateVariable         = subRegion.getField< rateAndState::stateVariable >();
        arrayView2d< real64 > const deltaSlip             = subRegion.getField< rateAndState::deltaSlip >();
        arrayView2d< real64 > const dispJump              = subRegion.getField< contact::dispJump >();
        arrayView2d< real64 > const rateStateError        = subRegion.getField< rateAndState::error >();

        // Stage rates
        arrayView2d< real64 const > const stateVariableRKStageRate  = subRegion.getField< rateAndState::stateVariableRKStageRate >();
        arrayView3d< real64 const > const deltaSlipRKStageRate      = subRegion.getField< rateAndState::deltaSlipRKStageRate >();

        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {   
            // Runge-Kutta stage rates for stage 3
            real64 const stateVariableRKStageRate3 = frictionKernelWrapper.stateEvolution(k, slipRate[k], stateVariable[k]);
            real64 const deltaSlipRKStageRate3[2]  = {slipVelocity[k][0], slipVelocity[k][1]};
            
            // Third order update
            stateVariable[k] = stateVariable_n[k] + dtAdaptive/6*( stateVariableRKStageRate[k][0] + 4*stateVariableRKStageRate[k][1] + stateVariableRKStageRate3 );
            deltaSlip[k][0]  = deltaSlip_n[k][0]  + dtAdaptive/6*( deltaSlipRKStageRate[k][0][0]  + 4*deltaSlipRKStageRate[k][1][0] + deltaSlipRKStageRate3[0] );
            deltaSlip[k][1]  = deltaSlip_n[k][1]  + dtAdaptive/6*( deltaSlipRKStageRate[k][0][1]  + 4*deltaSlipRKStageRate[k][1][1] + deltaSlipRKStageRate3[1] );
            
            dispJump[k][1] = dispJump_n[k][1] + deltaSlip[k][0];
            dispJump[k][2] = dispJump_n[k][2] + deltaSlip[k][1];

            // Second order updates used for adaptive error comparison
            real64 const stateVariableRK2 = stateVariable_n[k] + dtAdaptive/2*(stateVariableRKStageRate[k][0] + stateVariableRKStageRate3);
            real64 const deltaSlipRK2[2] = { deltaSlip_n[k][0] + dtAdaptive/2*( deltaSlipRKStageRate[k][0][1]  + deltaSlipRKStageRate3[0] ),
                                             deltaSlip_n[k][1] + dtAdaptive/2*( deltaSlipRKStageRate[k][0][1]  + deltaSlipRKStageRate3[1] ) }; 
            
            // TODO: Change error computation
            // Compute relative errors based on the maximum error in state or slip
            
            rateStateError[k][0] = (deltaSlip[k][0] - deltaSlipRK2[0]) / 
                                   ( m_timestepAbsTol + m_timestepRelTol * LvArray::math::max( LvArray::math::abs(deltaSlip[k][0]), LvArray::math::abs(deltaSlipRK2[0]) ));
            rateStateError[k][1] = (deltaSlip[k][1] - deltaSlipRK2[1]) / 
                                   ( m_timestepAbsTol + m_timestepRelTol * LvArray::math::max( LvArray::math::abs(deltaSlip[k][1]), LvArray::math::abs(deltaSlipRK2[1]) ));
            rateStateError[k][2] = (stateVariable[k] - stateVariableRK2) / 
                                   ( m_timestepAbsTol + m_timestepRelTol * LvArray::math::max( LvArray::math::abs(stateVariable[k]), LvArray::math::abs(stateVariableRK2) ));
            // real64 const errorStateVariable = LvArray::math::abs( (stateVariable[k] - stateVariableRK2) /( stateVariable[k] + std::numeric_limits< real64 >::epsilon() ) );
            // real64 const errorSlip[2] = { ( deltaSlip[k][0] - deltaSlipRK2[0] )/( deltaSlip[k][0] + std::numeric_limits< real64 >::epsilon() ),
            //                               ( deltaSlip[k][1] - deltaSlipRK2[1] )/( deltaSlip[k][0] + std::numeric_limits< real64 >::epsilon() ) };
        } );
      } );
    } );
    // Update timestep based on the time step error 
    dtAdaptive = setNextDt(dtAdaptive, domain);
    if (m_successfulStep) // set in setNextDt
    {
      // Successful step. Compute stresses, and slip velocity and save results, then exit the adaptive time step loop
      dtStress = updateStresses( time_n, dt, cycleNumber, domain );
      updateSlipVelocity( time_n, dt, domain );
      saveState(domain);
      break;
    }

  }
  
  // return time step size achieved by stress solver
  return dtAdaptive;
}

real64 QuasiDynamicEQRK32::updateStresses( real64 const & time_n,
                                            real64 const & dt,
                                            const int cycleNumber,
                                            DomainPartition & domain ) const
{
  GEOS_LOG_LEVEL_RANK_0( 1, "Stress solver" );
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
        arrayView2d< real64 const > const traction_n      = subRegion.getField< fields::contact::traction_n >();

        string const & fricitonLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
        RateAndStateFriction const & frictionLaw = getConstitutiveModel< RateAndStateFriction >( subRegion, fricitonLawName );

        RateAndStateFriction::KernelWrapper frictionKernelWrapper = frictionLaw.createKernelUpdates();

        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {
          SpringSliderParameters springSliderParameters = SpringSliderParameters( traction[k][0],
                                                                                  frictionKernelWrapper.getACoefficient( k ),
                                                                                  frictionKernelWrapper.getBCoefficient( k ),
                                                                                  frictionKernelWrapper.getDcCoefficient( k ) );


          traction[k][1] = traction_n[k][1] + springSliderParameters.tauRate * dt
                           - springSliderParameters.springStiffness * deltaSlip[k][0];
          traction[k][2] = traction_n[k][2] + springSliderParameters.tauRate * dt
                           - springSliderParameters.springStiffness * deltaSlip[k][1];
        } );
      } );
    } );
    return dt;
  }
}

void QuasiDynamicEQRK32::updateSlipVelocity( real64 const & time_n,
                                               real64 const & dt,
                                               DomainPartition & domain ) const
{
  GEOS_LOG_LEVEL_RANK_0( 1, "Rate and State solver" );
  integer const maxIterNewton = m_nonlinearSolverParameters.m_maxIterNewton;
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
      rateAndStateKernels::createAndLaunch<rateAndStateKernels::ExplicitRateAndStateKernel, parallelDevicePolicy<> >( subRegion, viewKeyStruct::frictionLawNameString(), m_shearImpedance, maxIterNewton, newtonTol, time_n, dt );
    } );
  } );
}

void QuasiDynamicEQRK32::saveState( DomainPartition & domain) const
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )

    {
      mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                             [&]( localIndex const,
                                                                                  SurfaceElementSubRegion & subRegion )
      {
        arrayView1d< real64 const > const stateVariable = subRegion.getField< rateAndState::stateVariable >();
        arrayView2d< real64 const > const slipVelocity  = subRegion.getField< rateAndState::slipVelocity >();
        arrayView2d< real64 const > const deltaSlip     = subRegion.getField< rateAndState::deltaSlip >();
        arrayView2d< real64 const > const dispJump      = subRegion.getField< contact::dispJump >();
        arrayView2d< real64 const > const traction      = subRegion.getField< contact::traction >();

        arrayView1d< real64 > const stateVariable_n = subRegion.getField< rateAndState::stateVariable_n >();
        arrayView2d< real64 > const slipVelocity_n  = subRegion.getField< rateAndState::slipVelocity_n >();
        arrayView2d< real64 > const deltaSlip_n     = subRegion.getField< rateAndState::deltaSlip >();
        arrayView2d< real64 > const dispJump_n      = subRegion.getField< contact::dispJump_n >();
        arrayView2d< real64 > const traction_n      = subRegion.getField< contact::traction_n >();

        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {
          stateVariable_n[k]  = stateVariable[k];
          LvArray::tensorOps::copy< 2 >( deltaSlip_n[k], deltaSlip[k] );
          LvArray::tensorOps::copy< 2 >( slipVelocity_n[k], slipVelocity[k] );
          LvArray::tensorOps::copy< 3 >( dispJump_n[k], dispJump[k] );
          LvArray::tensorOps::copy< 3 >( traction_n[k], traction[k] );
        } );
      } );
    } );
}

real64 QuasiDynamicEQRK32::setNextDt( real64 const & currentDt, DomainPartition & domain )
{

  // Spring-slider shear traction computation
  real64 scaledL2Error  = 0.0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion const & subRegion )
    {
      arrayView2d< real64 const > const error = subRegion.getField< rateAndState::error >();

      RAJA::ReduceSum< parallelDeviceReduce, real64 > scaledl2ErrorSquared( 0.0 );
      integer const N = subRegion.size();
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        scaledl2ErrorSquared += LvArray::tensorOps::l2NormSquared<3>(error[k]);
      } );
      scaledL2Error = LvArray::math::sqrt(MpiWrapper::sum( scaledl2ErrorSquared.get() / (3.0*N) ));
  } );
  } );

  // PID error controller + limiter
  real64 const k = LvArray::math::min(m_rkOrders[0], m_rkOrders[1]) + 1.0;
  real64 const eps0 = 1.0/(scaledL2Error + std::numeric_limits< real64 >::epsilon()); // n + 1
  real64 const eps1 =  1.0/(m_prevTimestepErrors[0] + std::numeric_limits< real64 >::epsilon()); // n
  real64 const eps2 = 1.0/(m_prevTimestepErrors[1] + std::numeric_limits< real64 >::epsilon()); // n-1
  // Limiter is 1.0 + atan(x - 1.0). Here use atan(x) = atan2(x, 1.0)
  real64 const dtFactor = 1.0 + LvArray::math::atan2( pow(eps0, m_beta[0] / k ) *  pow(eps1, m_beta[1] / k ) *  pow(eps2, m_beta[2] / k ) - 1.0, 1.0);
  
  real64 const nextDt = dtFactor*currentDt;
  m_successfulStep = (dtFactor >= m_timestepAcceptSafety) ? true : false;
  if ( m_successfulStep )
  {
    m_prevTimestepErrors[1] = m_prevTimestepErrors[0];
    m_prevTimestepErrors[0] = scaledL2Error;
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "Adaptive time step successful. The next dt will be {:.2e} s", nextDt ));  
  }
  else
  {
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "Adaptive time step failed. The next dt will be {:.2e} s", nextDt ));  
  }

  return nextDt;
}

REGISTER_CATALOG_ENTRY( SolverBase, QuasiDynamicEQRK32, string const &, dataRepository::Group * const )

} // namespace geos
