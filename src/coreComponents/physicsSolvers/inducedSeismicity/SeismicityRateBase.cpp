/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SeismicityRateBase.cpp
 */

#include "SeismicityRateBase.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

//START_SPHINX_INCLUDE_CONSTRUCTOR
SeismicityRateBase::SeismicityRateBase( const string & name,
                              Group * const parent ):
  SolverBase( name, parent ),
  m_stressSolver( nullptr )
  {
    this->registerWrapper( viewKeyStruct::stressSolverNameString(), &m_stressSolverName ).
          setInputFlag( InputFlags::REQUIRED ).
          setDescription( "Name of solver for computing stress" );
    this->registerWrapper( viewKeyStruct::initialFaultNormalStressString(), &m_initialFaultNormalStress ).
          setInputFlag( InputFlags::OPTIONAL ).
          setDescription( "Initial normal stress" );
    this->registerWrapper( viewKeyStruct::initialFaultShearStressString(), &m_initialFaultShearStress ).
          setInputFlag( InputFlags::OPTIONAL ).
          setDescription( "Initial shear stress" );

    // Only temporarily store user-specified fault directions here, then initialize fault projection tensors
    this->registerWrapper( viewKeyStruct::faultNormalString(), &m_faultNormalDirection ).
          setInputFlag( InputFlags::OPTIONAL ).
          setDescription( "Fault normal direction" );
    this->registerWrapper( viewKeyStruct::faultShearString(), &m_faultShearDirection ).
          setInputFlag( InputFlags::OPTIONAL ).
          setDescription( "Fault shear direction" );

    // Check orthogonality of user-specified faults
    if ( std::abs(LvArray::tensorOps::AiBi< 3 >(m_faultNormalDirection, m_faultShearDirection)) > 1e-8 )
    {
      GEOS_ERROR("Fault normal and fault shear must be orthogonal");
    }
    // If user has specified faults, normalize them to be unit vectors
    else if ( LvArray::tensorOps::l2Norm< 3 >( m_faultNormalDirection ) > 1e-8 )
    {
      LvArray::tensorOps::normalize< 3 >( m_faultNormalDirection );
      LvArray::tensorOps::normalize< 3 >( m_faultShearDirection );
    }

    
  }
//END_SPHINX_INCLUDE_CONSTRUCTOR

//START_SPHINX_INCLUDE_REGISTERDATAONMESH
void SeismicityRateBase::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< inducedSeismicity::initialMeanNormalStress >( getName() );
      subRegion.registerField< inducedSeismicity::initialMeanShearStress >( getName() );

      subRegion.registerField< inducedSeismicity::meanNormalStress >( getName() );
      subRegion.registerField< inducedSeismicity::meanNormalStress_n >( getName() );
      subRegion.registerField< inducedSeismicity::meanShearStress >( getName() );
      subRegion.registerField< inducedSeismicity::meanShearStress_n >( getName() );

      subRegion.registerField< inducedSeismicity::seismicityRate >( getName() );
    } );
   } );
}

void SeismicityRateBase::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const tempR = subRegion.getField< inducedSeismicity::seismicityRate >();
      tempR.setValues< parallelHostPolicy >( 1.0 );
    } );
  } );
}

void SeismicityRateBase::updateMeanSolidStress( ElementSubRegionBase & subRegion )
{
  // Retrieve field variables
  arrayView1d< real64 > const sig = subRegion.getField< inducedSeismicity::meanNormalStress >();
  arrayView1d< real64 > const sig_n = subRegion.getField< inducedSeismicity::meanNormalStress_n >();
  arrayView1d< real64 > const tau = subRegion.getField< inducedSeismicity::meanShearStress >();
  arrayView1d< real64 > const tau_n = subRegion.getField< inducedSeismicity::meanShearStress_n >();

  // Retrieve stress state computed by m_stressSolver, called in solverStep
  string const & solidModelName = subRegion.getReference< string >( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString());
  SolidBase const & solidModel = getConstitutiveModel< SolidBase >( subRegion, solidModelName );
  arrayView3d< real64 const, solid::STRESS_USD > const stress = solidModel.getStress();

  // Construct Voigt notation projection tensor (dyadic product of fault normal and shear vectors)
  // for when multiplying in inner dot product with symmetric stress tensor
  real64 faultNormalProjectionTensor[6];
  real64 faultShearProjectionTensor[6];
  constructFaultStressProjectionTensors(faultNormalProjectionTensor, faultShearProjectionTensor);

  // loop over all elements and computed unweighted average across all nodes 
  // as average stress state acting over the cell ceneter
  // TODO: APPLY WEIGHTS TO AVERAGE TO ACCOUNT FOR ELEMENT SHAPE
  forAll< parallelDevicePolicy<> >(  sig.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    // update projected stresses from previou step
    sig_n[k] = sig[k];
    tau_n[k] = tau[k];

    // Compute average stress state
    array1d< real64 > meanStress( 6 );
    for( int q = 0; q < stress.size(1); q++ ) 
    {
      LvArray::tensorOps::add< 6 >( meanStress, stress[k][q] );
    }
    LvArray::tensorOps::scale< 6 >(meanStress, 1./stress.size(1));

    // Project average stress state to fault
    sig[k] = LvArray::tensorOps::AiBi< 6 >( meanStress, faultNormalProjectionTensor);
    tau[k] = LvArray::tensorOps::AiBi< 6 >( meanStress, faultShearProjectionTensor);
  } );
}

void SeismicityRateBase::initializeMeanSolidStress(real64 const time_n, integer const cycleNumber, DomainPartition & domain)
{
  // Only call initialization step before stress solver has been called for first time step
  if ( cycleNumber == 0 ) {
    // Call solverStep of stress solver with dt=0 to initialize stresses in the matrix
    m_stressSolver->solverStep(time_n, 0.0, cycleNumber, domain );

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )

    {
      mesh.getElemManager().forElementSubRegions( regionNames,
                                                  [&]( localIndex const,
                                                       ElementSubRegionBase & subRegion )
      {
        // Retrieve field variables
        arrayView1d< real64 const > const sig = subRegion.getField< inducedSeismicity::meanNormalStress >();
        arrayView1d< real64 > const sig_i = subRegion.getField< inducedSeismicity::initialMeanNormalStress >();
        arrayView1d< real64 const > const tau = subRegion.getField< inducedSeismicity::meanShearStress >();
        arrayView1d< real64 > const tau_i = subRegion.getField< inducedSeismicity::initialMeanShearStress >();

        // If solid stress solver has been called, call updateMeanSolidStress to project stress state onto faults
        if ( subRegion.hasWrapper( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString() ) )
        {
          if ( LvArray::tensorOps::l2Norm< 3 >( m_faultNormalDirection ) < 1e-8 )
          {
            GEOS_ERROR("Fault directions must be specified for solid stress solver");
          }

          updateMeanSolidStress( subRegion );
          forAll< parallelDevicePolicy<> >(  sig.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
          {
            // Set initial stress conditions on faults
            sig_i[k] = sig[k];
            tau_i[k] = tau[k];
          } );

        }
        // If solid stress solver has not been called, initialize fault stresses to user-defined values
        else
        {
          arrayView1d< real64 > const tempSigIni = subRegion.getField< inducedSeismicity::initialMeanNormalStress >();
          tempSigIni.setValues< parallelHostPolicy >( m_initialFaultNormalStress );
          arrayView1d< real64 > const tempSig = subRegion.getField< inducedSeismicity::meanNormalStress >();
          tempSig.setValues< parallelHostPolicy >( m_initialFaultNormalStress );
          arrayView1d< real64 > const tempSig_n = subRegion.getField< inducedSeismicity::meanNormalStress_n >();
          tempSig_n.setValues< parallelHostPolicy >( m_initialFaultNormalStress );

          arrayView1d< real64 > const tempTauIni = subRegion.getField< inducedSeismicity::initialMeanShearStress >();
          tempTauIni.setValues< parallelHostPolicy >( m_initialFaultShearStress );
          arrayView1d< real64 > const tempTau = subRegion.getField< inducedSeismicity::meanShearStress >();
          tempTau.setValues< parallelHostPolicy >( m_initialFaultShearStress );
          arrayView1d< real64 > const tempTau_n = subRegion.getField< inducedSeismicity::meanShearStress_n >();
          tempTau_n.setValues< parallelHostPolicy >( m_initialFaultShearStress );
        }
      } );
    } );
  }
}

void SeismicityRateBase::constructFaultStressProjectionTensors( 
  real64 faultNormalProjectionTensor[6], real64 faultShearProjectionTensor[6] )
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

void SeismicityRateBase::postProcessInput()
{
  // Initialize member stress solver as specified in XML input
  m_stressSolver = &this->getParent().getGroup< SolverBase >( m_stressSolverName );

  SolverBase::postProcessInput();
}

SeismicityRateBase::~SeismicityRateBase()
{
  // TODO Auto-generated destructor stub
}

} // namespace geos
