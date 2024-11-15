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
 * @file SeismicityRate.cpp
 */

#include "SeismicityRate.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "kernels/SeismicityRateKernels.hpp"
#include "physicsSolvers/inducedSeismicity/inducedSeismicityFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

#include "fieldSpecification/FieldSpecificationManager.hpp"


namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

SeismicityRate::SeismicityRate( const string & name,
                                Group * const parent ):
  PhysicsSolverBase( name, parent ),
  m_stressSolver( nullptr )
{
  this->registerWrapper( viewKeyStruct::directEffectString(), &m_directEffect ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Rate-and-state friction direct effect parameter." );

  this->registerWrapper( viewKeyStruct::backgroundStressingRateString(), &m_backgroundStressingRate ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Background stressing rate (Pa/s)." );

  this->registerWrapper( viewKeyStruct::stressSolverNameString(), &m_stressSolverName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of solver for computing stress" );

  this->registerWrapper( viewKeyStruct::faultNormalDirectionString(), &m_faultNormalDirection ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fault normal direction" );

  this->registerWrapper( viewKeyStruct::faultShearDirectionString(), &m_faultShearDirection ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fault shear direction" );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );
}

void SeismicityRate::postInputInitialization()
{
  // Check orthogonality of user-specified faults
  if( std::abs( LvArray::tensorOps::AiBi< 3 >( m_faultNormalDirection, m_faultShearDirection )) > 1e-8 )
  {
    GEOS_ERROR( "Fault normal and fault shear directions must be orthogonal" );
  }
  // If user has specified faults, normalize them to be unit vectors
  else if( LvArray::tensorOps::l2Norm< 3 >( m_faultNormalDirection ) > 1e-8 )
  {
    LvArray::tensorOps::normalize< 3 >( m_faultNormalDirection );
    LvArray::tensorOps::normalize< 3 >( m_faultShearDirection );
  }

  // Initialize member stress solver as specified in XML input
  if( !m_stressSolverName.empty() )
  {
    m_stressSolver = &this->getParent().getGroup< PhysicsSolverBase >( m_stressSolverName );
  }

  PhysicsSolverBase::postInputInitialization();
}

SeismicityRate::~SeismicityRate()
{
  // TODO Auto-generated destructor stub
}

void SeismicityRate::registerDataOnMesh( Group & meshBodies )
{
  PhysicsSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< inducedSeismicity::initialProjectedNormalTraction >( getName() );

      subRegion.registerField< inducedSeismicity::initialProjectedShearTraction >( getName() );

      subRegion.registerField< inducedSeismicity::projectedNormalTraction >( getName() );

      subRegion.registerField< inducedSeismicity::projectedNormalTraction_n >( getName() );

      subRegion.registerField< inducedSeismicity::projectedShearTraction >( getName() );

      subRegion.registerField< inducedSeismicity::projectedShearTraction_n >( getName() );

      subRegion.registerField< inducedSeismicity::seismicityRate >( getName() );

      subRegion.registerField< inducedSeismicity::logDenom >( getName() );
    } );
  } );
}

void SeismicityRate::updateFaultTraction( ElementSubRegionBase & subRegion ) const
{
  // Retrieve field variables
  arrayView1d< real64 > const sig   = subRegion.getField< inducedSeismicity::projectedNormalTraction >();
  arrayView1d< real64 > const tau   = subRegion.getField< inducedSeismicity::projectedShearTraction >();

  // Retrieve stress state computed by m_stressSolver, called in solverStep
  string const & solidModelName = subRegion.getReference< string >( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString());
  SolidBase const & solidModel  = getConstitutiveModel< SolidBase >( subRegion, solidModelName );
  arrayView3d< real64 const, solid::STRESS_USD > const stress = solidModel.getStress();

  // Construct Voigt notation projection tensor (dyadic product of fault normal and shear vectors)
  // for when multiplying in inner dot product with symmetric stress tensor
  real64 faultNormalProjectionTensor[6]{};
  real64 faultShearProjectionTensor[6]{};
  constructFaultStressProjectionTensors( faultNormalProjectionTensor, faultShearProjectionTensor );

  // loop over all elements and compute unweighted average across all nodes
  // as average stress state acting over the cell ceneter
  // TODO: APPLY WEIGHTS TO AVERAGE TO ACCOUNT FOR ELEMENT SHAPE
  forAll< parallelDevicePolicy<> >( sig.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    // Compute average stress state
    real64 meanStress[ 6 ]{};
    for( int q = 0; q < stress.size( 1 ); q++ )
    {
      LvArray::tensorOps::add< 6 >( meanStress, stress[k][q] );
    }
    LvArray::tensorOps::scale< 6 >( meanStress, 1./stress.size( 1 ));

    // Project average stress state to fault
    sig[k] = LvArray::tensorOps::AiBi< 6 >( meanStress, faultNormalProjectionTensor );
    tau[k] = LvArray::tensorOps::AiBi< 6 >( meanStress, faultShearProjectionTensor );
  } );

  // For poroelastic models, we must calculate the total stress before computing the effective stresses
  // on the fault. This requires retrieving both the pressure field and the Biot coefficient. We first check
  // to see if a flow solver exists, retrive the pressure field, then pass the porous model through the lambda
  // cast to access the Biot coefficient. Finally, effective stresses on the fault are calculated.
  if( subRegion.hasWrapper( FlowSolverBase::viewKeyStruct::fluidNamesString() ) )
  {
    arrayView1d< real64 const > const pres = subRegion.getField< flow::pressure >();

    string const & porousSolidModelName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
    CoupledSolidBase & porousSolid = getConstitutiveModel< CoupledSolidBase >( subRegion, porousSolidModelName );
    constitutive::ConstitutivePassThru< CoupledSolidBase >::execute( porousSolid, [&] ( auto & castedPorousSolid )
    {
      // Initialize biotCoefficient as const arrayView before passing it through the lambda cast
      arrayView1d< real64 const > const biotCoefficient = castedPorousSolid.getBiotCoefficient();

      // To calculate the action of the total stress on the fault from our previous calculations,
      // we need to project the action of the pore pressure on the stress tensor onto the fault
      computeTotalStressOnFault( biotCoefficient,
                                 pres,
                                 faultNormalProjectionTensor,
                                 faultShearProjectionTensor,
                                 sig,
                                 tau );
    } );
  }
}

void SeismicityRate::computeTotalStressOnFault( arrayView1d< real64 const > const & biotCoefficient,
                                                arrayView1d< real64 const > const & pres,
                                                real64 const (&faultNormalProjectionTensor)[6],
                                                real64 const (&faultShearProjectionTensor)[6],
                                                arrayView1d< real64 > const & sig,
                                                arrayView1d< real64 > const & tau ) const
{
  // To calculate the action of the total stress on the fault from our previous calculations,
  // we need to project the action of the pore pressure on the stress tensor onto the fault
  forAll< parallelDevicePolicy<> >( sig.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    // Form pressure as tensor
    real64 pressureTensor[ 6 ]{};
    LvArray::tensorOps::symAddIdentity< 3 >( pressureTensor, -biotCoefficient[k]*pres[k] );

    // Project pressure tensor onto fault orientations
    real64 const pressureOnFaultNormal = LvArray::tensorOps::AiBi< 6 >( pressureTensor, faultNormalProjectionTensor );
    real64 const pressureOnFaultShear  = LvArray::tensorOps::AiBi< 6 >( pressureTensor, faultShearProjectionTensor );

    // Calculate total stress on the faults
    sig[k] += pressureOnFaultNormal;
    tau[k] += pressureOnFaultShear;
  } );

}

void SeismicityRate::initializeFaultTraction( real64 const time_n, integer const cycleNumber, DomainPartition & domain ) const
{
  // Only call initialization step before stress solver has been called for first time step
  if( cycleNumber == 0 )
  {
    // Call solverStep of stress solver with dt=0 to initialize stresses in the matrix

    updateStresses( time_n, 0.0, cycleNumber, domain );

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )

    {
      mesh.getElemManager().forElementSubRegions( regionNames,
                                                  [&]( localIndex const,
                                                       ElementSubRegionBase & subRegion )
      {

        saveOldState( subRegion );

        // Retrieve field variables
        arrayView1d< real64 const > const sig = subRegion.getField< inducedSeismicity::projectedNormalTraction >();
        arrayView1d< real64 const > const tau = subRegion.getField< inducedSeismicity::projectedShearTraction >();

        arrayView1d< real64 > const sig_i = subRegion.getField< inducedSeismicity::initialProjectedNormalTraction >();
        arrayView1d< real64 > const tau_i = subRegion.getField< inducedSeismicity::initialProjectedShearTraction >();

        forAll< parallelDevicePolicy<> >( sig.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {
          // Set initial stress conditions on faults
          sig_i[k] = sig[k];
          tau_i[k] = tau[k];
        } );
      } );
    } );
  }
}

void SeismicityRate::constructFaultStressProjectionTensors( real64 (& faultNormalProjectionTensor)[6],
                                                            real64 (& faultShearProjectionTensor)[6] ) const
{
  faultNormalProjectionTensor[0] = m_faultNormalDirection[0]*m_faultNormalDirection[0];
  faultNormalProjectionTensor[1] = m_faultNormalDirection[1]*m_faultNormalDirection[1];
  faultNormalProjectionTensor[2] = m_faultNormalDirection[2]*m_faultNormalDirection[2];
  faultNormalProjectionTensor[3] = 2*m_faultNormalDirection[1]*m_faultNormalDirection[2];
  faultNormalProjectionTensor[4] = 2*m_faultNormalDirection[0]*m_faultNormalDirection[2];
  faultNormalProjectionTensor[5] = 2*m_faultNormalDirection[0]*m_faultNormalDirection[1];

  faultShearProjectionTensor[0] = m_faultShearDirection[0]*m_faultNormalDirection[0];
  faultShearProjectionTensor[1] = m_faultShearDirection[1]*m_faultNormalDirection[1];
  faultShearProjectionTensor[2] = m_faultShearDirection[2]*m_faultNormalDirection[2];
  faultShearProjectionTensor[3] = m_faultShearDirection[1]*m_faultNormalDirection[2] + m_faultShearDirection[2]*m_faultNormalDirection[1];
  faultShearProjectionTensor[4] = m_faultShearDirection[0]*m_faultNormalDirection[2] + m_faultShearDirection[2]*m_faultNormalDirection[0];
  faultShearProjectionTensor[5] = m_faultShearDirection[0]*m_faultNormalDirection[1] + m_faultShearDirection[1]*m_faultNormalDirection[0];
}

real64 SeismicityRate::solverStep( real64 const & time_n,
                                   real64 const & dt,
                                   const int cycleNumber,
                                   DomainPartition & domain )
{
  // Save initial stress state on pre-defined fault orientations to field variables
  if( cycleNumber == 0 )
  {
    initializeFaultTraction( time_n, cycleNumber, domain );
  }

  real64 const dtStress = updateStresses( time_n, dt, cycleNumber, domain );

  // Loop over subRegions to solve for seismicity rate
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      // solve for the seismicity rate given new stresses on faults
      integralSolverStep( time_n, dtStress, subRegion );

      // save old state
      saveOldState( subRegion );
    } );
  } );

  // return time step size achieved by stress solver
  return dtStress;
}

real64 SeismicityRate::updateStresses( real64 const & time_n,
                                       real64 const & dt,
                                       const int cycleNumber,
                                       DomainPartition & domain ) const
{
  // Call member variable stress solver to update the stress state
  if( m_stressSolver )
  {

    // 1. Solve the momentum balance
    real64 const dtStress =  m_stressSolver->solverStep( time_n, dt, cycleNumber, domain );

    // 2. Loop over subRegions to update stress on faults
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )

    {
      mesh.getElemManager().forElementSubRegions( regionNames,
                                                  [&]( localIndex const,
                                                       ElementSubRegionBase & subRegion )
      {
        updateFaultTraction( subRegion );
      } );
    } );
    return dtStress;
  }
  else
  {
    char const bcLogMessage[] =
      "SeismicityRate {}: at time {}s, "
      "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. "
      "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
      "\nThe total number of target elements (including ghost elements) is {}. "
      "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set.";

    forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&] ( string const &,
                                                                   MeshLevel & mesh,
                                                                   arrayView1d< string const > const & )
    {

      FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

      std::vector< string > const keys = { inducedSeismicity::projectedNormalTraction::key(),
                                           inducedSeismicity::projectedShearTraction::key() };

      for( auto const & key : keys )
      {
        fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                                 mesh,
                                                 key,
                                                 [&]( FieldSpecificationBase const & fs,
                                                      string const & setName,
                                                      SortedArrayView< localIndex const > const & lset,
                                                      ElementSubRegionBase & subRegion,
                                                      string const & )
        {
          if( fs.getLogLevel() >= 1 )
          {
            globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( lset.size() );
            GEOS_LOG_RANK_0( GEOS_FMT( bcLogMessage,
                                       this->getName(), time_n+dt, FieldSpecificationBase::catalogName(),
                                       fs.getName(), setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
          }

          // Specify the bc value of the field
          fs.applyFieldValue< FieldSpecificationEqual,
                              parallelDevicePolicy<> >( lset,
                                                        time_n + dt,
                                                        subRegion,
                                                        key );
        } );
      }
    } );
  }
  return dt;
}

void SeismicityRate::saveOldState( ElementSubRegionBase & subRegion ) const
{
  // Retrieve field variables
  arrayView1d< real64 > const sig   = subRegion.getField< inducedSeismicity::projectedNormalTraction >();
  arrayView1d< real64 > const sig_n = subRegion.getField< inducedSeismicity::projectedNormalTraction_n >();
  arrayView1d< real64 > const tau   = subRegion.getField< inducedSeismicity::projectedShearTraction >();
  arrayView1d< real64 > const tau_n = subRegion.getField< inducedSeismicity::projectedShearTraction_n >();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    // update projected stresses from previou step
    sig_n[k] = sig[k];
    tau_n[k] = tau[k];
  } );
}

// Solve integral solution to ODE
void SeismicityRate::integralSolverStep( real64 const & time_n,
                                         real64 const & dt,
                                         ElementSubRegionBase & subRegion )
{
  if( subRegion.hasWrapper( FlowSolverBase::viewKeyStruct::fluidNamesString() )  )
  {
    seismicityRateKernels::createAndLaunch< parallelDevicePolicy<>, true >( subRegion, time_n, dt, m_directEffect, m_backgroundStressingRate );
  }
  else
  {
    seismicityRateKernels::createAndLaunch< parallelDevicePolicy<>, false >( subRegion, time_n, dt, m_directEffect, m_backgroundStressingRate );
  }
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SeismicityRate, string const &, dataRepository::Group * const )
} // namespace geos
