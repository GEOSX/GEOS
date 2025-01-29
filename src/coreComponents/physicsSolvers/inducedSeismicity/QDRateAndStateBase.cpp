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
 * @file QDRateAndStateBase.cpp
 */

#include "QDRateAndStateBase.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mesh/DomainPartition.hpp"
#include "rateAndStateFields.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "constitutive/contact/RateAndStateFriction.hpp"
#include "kernels/RateAndStateKernelsBase.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

QDRateAndStateBase::QDRateAndStateBase( const string & name,
                                        Group * const parent ):
  PhysicsSolverBase( name, parent ),
  m_shearImpedance( 0.0 )
{
  this->registerWrapper( viewKeyStruct::shearImpedanceString(), &m_shearImpedance ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear impedance." );
}

QDRateAndStateBase::~QDRateAndStateBase()
{
  // TODO Auto-generated destructor stub
}

void QDRateAndStateBase::registerDataOnMesh( Group & meshBodies )
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
      subRegion.registerField< rateAndState::slipVelocity_n >( getName() ).
        setDimLabels( 1, labels2Comp ).reference().resizeDimension< 1 >( 2 );

      subRegion.registerField< rateAndState::shearTraction >( getName() ).
        setDimLabels( 1, labels2Comp ).reference().resizeDimension< 1 >( 2 );
      subRegion.registerField< rateAndState::normalTraction >( getName() );

      subRegion.registerField< rateAndState::shearTraction_n >( getName() ).
        setDimLabels( 1, labels2Comp ).reference().resizeDimension< 1 >( 2 );
      subRegion.registerField< rateAndState::normalTraction_n >( getName() );

      subRegion.registerField< rateAndState::backgroundShearStress >( getName() ).
        setDimLabels( 1, labels2Comp ).reference().resizeDimension< 1 >( 2 );
      subRegion.registerField< rateAndState::backgroundNormalStress >( getName() );
    } );
  } );
}

void QDRateAndStateBase::enforceRateAndVelocityConsistency( SurfaceElementSubRegion & subRegion ) const
{
  arrayView2d< real64 > const slipVelocity = subRegion.getField< rateAndState::slipVelocity >();
  arrayView1d< real64 > const slipRate  = subRegion.getField< rateAndState::slipRate >();
  arrayView1d< real64 const > const stateVariable  = subRegion.getField< rateAndState::stateVariable >();

  arrayView2d< real64 > const backgroundShearStress = subRegion.getField< rateAndState::backgroundShearStress >();
  arrayView1d< real64 > const backgroundNormalStress = subRegion.getField< rateAndState::backgroundNormalStress >();

  real64 const shearImpedance = m_shearImpedance;

  string const & frictionaLawName = subRegion.getReference< string >( viewKeyStruct::frictionLawNameString() );
  constitutive::RateAndStateFriction const & frictionLaw = subRegion.getConstitutiveModel< constitutive::RateAndStateFriction >( frictionaLawName );

  constitutive::RateAndStateFriction::KernelWrapper frictionLawKernelWrapper = frictionLaw.createKernelUpdates();

  RAJA::ReduceMax< parallelDeviceReduce, int > negativeSlipRate( 0 );
  RAJA::ReduceMax< parallelDeviceReduce, int > bothNonZero( 0 );

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    if( slipRate[k] < 0.0 )
    {
      negativeSlipRate.max( 1 );
    }
    else if( LvArray::tensorOps::l2Norm< 2 >( slipVelocity[k] ) > 0.0 && slipRate[k] > 0.0 )
    {
      bothNonZero.max( 1 );
    }
    else if( LvArray::tensorOps::l2Norm< 2 >( slipVelocity[k] ) > 0.0 )
    {
      slipRate[k] = LvArray::tensorOps::l2Norm< 2 >( slipVelocity[k] );
    }
    else if( slipRate[k] > 0.0 )
    {
      real64 const frictionCoefficient = frictionLawKernelWrapper.frictionCoefficient( k, slipRate[k], stateVariable[k] );
      rateAndStateKernels::projectSlipRateBase( k,
                                                frictionCoefficient,
                                                shearImpedance,
                                                backgroundNormalStress,
                                                backgroundShearStress,
                                                slipRate,
                                                slipVelocity );
    }
  } );


  GEOS_ERROR_IF( negativeSlipRate.get() > 0, "SlipRate cannot be negative." );
  GEOS_ERROR_IF( bothNonZero.get() > 0, "Only one between slipRate and slipVelocity can be specified as i.c." );

}

void QDRateAndStateBase::applyInitialConditionsToFault( int const cycleNumber,
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
        enforceRateAndVelocityConsistency( subRegion );

        arrayView1d< real64 const > const slipRate        = subRegion.getField< rateAndState::slipRate >();
        arrayView1d< real64 const > const stateVariable   = subRegion.getField< rateAndState::stateVariable >();

        arrayView1d< real64 > const normalTraction  = subRegion.getField< rateAndState::normalTraction >();
        arrayView2d< real64 > const shearTraction   = subRegion.getField< rateAndState::shearTraction >();

        arrayView1d< real64 > const stateVariable_n       = subRegion.getField< rateAndState::stateVariable_n >();
        arrayView1d< real64 > const slipRate_n            = subRegion.getField< rateAndState::slipRate_n >();
        arrayView1d< real64 > const normalTraction_n      = subRegion.getField< rateAndState::normalTraction_n >();
        arrayView2d< real64 > const shearTraction_n       = subRegion.getField< rateAndState::shearTraction_n >();

        arrayView2d< real64 const > const backgroundShearStress  = subRegion.getField< rateAndState::backgroundShearStress >();
        arrayView1d< real64 const > const backgroundNormalStress = subRegion.getField< rateAndState::backgroundNormalStress >();

        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {
          slipRate_n [k]       = slipRate[k];
          stateVariable_n[k]   = stateVariable[k];

          normalTraction[k] = backgroundNormalStress[k];
          LvArray::tensorOps::copy< 2 >( shearTraction[k], backgroundShearStress[k] );

          normalTraction_n[k]  = normalTraction[k];
          LvArray::tensorOps::copy< 2 >( shearTraction_n[k], shearTraction[k] );
        } );
      } );
    } );
  }
}

void QDRateAndStateBase::saveState( DomainPartition & domain ) const
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
      arrayView2d< real64 const > const deltaSlip     = subRegion.getField< contact::deltaSlip >();
      arrayView2d< real64 const > const dispJump      = subRegion.getField< contact::dispJump >();
      arrayView2d< real64 const > const shearTraction      = subRegion.getField< rateAndState::shearTraction >();
      arrayView1d< real64 const > const normalTraction      = subRegion.getField< rateAndState::normalTraction >();


      arrayView1d< real64 > const stateVariable_n = subRegion.getField< rateAndState::stateVariable_n >();
      arrayView2d< real64 > const slipVelocity_n  = subRegion.getField< rateAndState::slipVelocity_n >();
      arrayView2d< real64 > const deltaSlip_n     = subRegion.getField< contact::deltaSlip >();
      arrayView2d< real64 > const dispJump_n      = subRegion.getField< contact::dispJump_n >();
      arrayView2d< real64 > const shearTraction_n = subRegion.getField< rateAndState::shearTraction_n >();
      arrayView1d< real64 > const normalTraction_n = subRegion.getField< rateAndState::normalTraction_n >();


      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        stateVariable_n[k]  = stateVariable[k];
        normalTraction_n[k] = normalTraction[k];
        LvArray::tensorOps::copy< 2 >( deltaSlip_n[k], deltaSlip[k] );
        LvArray::tensorOps::copy< 2 >( slipVelocity_n[k], slipVelocity[k] );
        LvArray::tensorOps::copy< 3 >( dispJump_n[k], dispJump[k] );
        LvArray::tensorOps::copy< 2 >( shearTraction_n[k], shearTraction[k] );
      } );
    } );
  } );
}

} // namespace geos
