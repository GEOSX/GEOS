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

/*
 * ContactSolverBase.cpp
 */

#include "ContactSolverBase.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/contact/FrictionBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "physicsSolvers/contact/LogLevelsInfo.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields::contact;

ContactSolverBase::ContactSolverBase( const string & name,
                                      Group * const parent ):
  SolidMechanicsLagrangianFEM( name, parent )
{
  this->getWrapper< string >( viewKeyStruct::contactRelationNameString() ).
    setInputFlag( dataRepository::InputFlags::FALSE );

  this->getWrapper< string >( viewKeyStruct::surfaceGeneratorNameString() ).
    setInputFlag( dataRepository::InputFlags::FALSE );
}

void ContactSolverBase::postInputInitialization()
{
  SolidMechanicsLagrangianFEM::postInputInitialization();

  GEOS_THROW_IF( m_timeIntegrationOption != TimeIntegrationOption::QuasiStatic,
                 GEOS_FMT( "{} {}: The attribute `{}` must be `{}`",
                           this->getCatalogName(), this->getName(),
                           viewKeyStruct::timeIntegrationOptionString(),
                           EnumStrings< TimeIntegrationOption >::toString( TimeIntegrationOption::QuasiStatic ) ),
                 InputError );
}

void ContactSolverBase::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  SolidMechanicsLagrangianFEM::registerDataOnMesh( meshBodies );

  setFractureRegions( meshBodies );

  forFractureRegionOnMeshTargets( meshBodies, [&] ( SurfaceElementRegion & fractureRegion )
  {
    string const labels[3] = { "normal", "tangent1", "tangent2" };

    fractureRegion.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion & subRegion )
    {
      setConstitutiveNamesCallSuper( subRegion );

      subRegion.registerField< fields::contact::dispJump >( getName() ).
        setDimLabels( 1, labels ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< fields::contact::deltaDispJump >( getName() ).
        setDimLabels( 1, labels ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< fields::contact::oldDispJump >( getName() ).
        setDimLabels( 1, labels ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< fields::contact::traction >( getName() ).
        setDimLabels( 1, labels ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< fields::contact::fractureState >( getName() );

      subRegion.registerField< fields::contact::oldFractureState >( getName() );

      subRegion.registerField< fields::contact::slip >( getName() );
    } );

  } );
}

void ContactSolverBase::setFractureRegions( dataRepository::Group const & meshBodies )
{
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel const & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementRegions< SurfaceElementRegion >( regionNames, [&] ( localIndex const, SurfaceElementRegion const & region )
    {
      m_fractureRegionNames.push_back( region.getName() );
    } );
  } );

  // TODO remove once multiple regions are fully supported
  GEOS_THROW_IF( m_fractureRegionNames.size() > 1,
                 GEOS_FMT( "{} {}: The number of fracture regions can not be more than one",
                           this->getCatalogName(), this->getName() ),
                 InputError );
}

void ContactSolverBase::computeFractureStateStatistics( MeshLevel const & mesh,
                                                        globalIndex & numStick,
                                                        globalIndex & numNewSlip,
                                                        globalIndex & numSlip,
                                                        globalIndex & numOpen ) const
{
  ElementRegionManager const & elemManager = mesh.getElemManager();

  array1d< globalIndex > localCounter( 4 );

  elemManager.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion const & subRegion )
  {
    arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
    arrayView1d< integer const > const & fractureState = subRegion.getField< fields::contact::fractureState >();

    RAJA::ReduceSum< parallelHostReduce, localIndex > stickCount( 0 ), newSlipCount( 0 ), slipCount( 0 ), openCount( 0 );
    forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
    {
      if( ghostRank[kfe] < 0 )
      {
        switch( fractureState[kfe] )
        {
          case FractureState::Stick:
            {
              stickCount += 1;
              break;
            }
          case FractureState::NewSlip:
            {
              newSlipCount += 1;
              break;
            }
          case FractureState::Slip:
            {
              slipCount += 1;
              break;
            }
          case FractureState::Open:
            {
              openCount += 1;
              break;
            }
        }
      }
    } );

    localCounter[0] += stickCount.get();
    localCounter[1] += newSlipCount.get();
    localCounter[2] += slipCount.get();
    localCounter[3] += openCount.get();
  } );

  array1d< globalIndex > totalCounter( 4 );

  MpiWrapper::allReduce( localCounter.data(),
                         totalCounter.data(),
                         4,
                         MPI_SUM,
                         MPI_COMM_GEOS );

  numStick    = totalCounter[0];
  numNewSlip  = totalCounter[1];
  numSlip     = totalCounter[2];
  numOpen     = totalCounter[3];
}

void ContactSolverBase::outputConfigurationStatistics( DomainPartition const & domain ) const
{
  if( getLogLevel() >=1 )
  {
    globalIndex numStick = 0;
    globalIndex numNewSlip  = 0;
    globalIndex numSlip  = 0;
    globalIndex numOpen  = 0;

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel const & mesh,
                                                                 arrayView1d< string const > const & )
    {
      computeFractureStateStatistics( mesh, numStick, numNewSlip, numSlip, numOpen );

      GEOS_LOG_RANK_0( GEOS_FMT( "  Number of element for each fracture state:"
                                 " stick: {:12} | new slip: {:12} | slip:  {:12} | open:  {:12}",
                                 numStick, numNewSlip, numSlip, numOpen ) );
    } );
  }
}

real64 ContactSolverBase::explicitStep( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                        real64 const & dt,
                                        const int GEOS_UNUSED_PARAM( cycleNumber ),
                                        DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  GEOS_MARK_FUNCTION;
  GEOS_ERROR( getDataContext() << ": ExplicitStep non available for contact solvers." );
  return dt;
}

void ContactSolverBase::synchronizeFractureState( DomainPartition & domain ) const
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { fields::contact::fractureState::key() }, { getUniqueFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors() );
  } );
}

void ContactSolverBase::setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const
{
  if( dynamic_cast< CellElementSubRegion * >( &subRegion ) )
  {
    SolidMechanicsLagrangianFEM::setConstitutiveNamesCallSuper( subRegion );
  }
  else if( dynamic_cast< SurfaceElementSubRegion * >( &subRegion ) )
  {
    subRegion.registerWrapper< string >( viewKeyStruct::frictionLawNameString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRestartFlags( RestartFlags::NO_WRITE ).
      setSizedFromParent( 0 );

    string & frictionLawName = subRegion.getReference< string >( viewKeyStruct::frictionLawNameString() );
    frictionLawName = PhysicsSolverBase::getConstitutiveName< FrictionBase >( subRegion );
    GEOS_ERROR_IF( frictionLawName.empty(), GEOS_FMT( "{}: FrictionBase model not found on subregion {}",
                                                      getDataContext(), subRegion.getDataContext() ) );
  }
}

} /* namespace geos */
