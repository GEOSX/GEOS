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
 * @file WaveSolverBase.cpp
 */

#include "WaveSolverBase.hpp"

#include "dataRepository/KeyNames.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"

#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/PerfectlyMatchedLayer.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/DomainPartition.hpp"
#include "WaveSolverUtils.hpp"
#include "events/EventManager.hpp"

#include <limits>

namespace geos
{

using namespace dataRepository;

WaveSolverBase::WaveSolverBase( const std::string & name,
                                Group * const parent ):
  SolverBase( name,
              parent )
{

  registerWrapper( viewKeyStruct::sourceCoordinatesString(), &m_sourceCoordinates ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Coordinates (x,y,z) of the sources" );

  registerWrapper( viewKeyStruct::receiverCoordinatesString(), &m_receiverCoordinates ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Coordinates (x,y,z) of the receivers" );

  registerWrapper( viewKeyStruct::timeSourceDelayString(), &m_timeSourceDelay ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1 ).
    setDescription( "Source time delay (1 / f0 by default)" );

  registerWrapper( viewKeyStruct::timeSourceFrequencyString(), &m_timeSourceFrequency ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Central frequency for the time source" );

  registerWrapper( viewKeyStruct::rickerOrderString(), &m_rickerOrder ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 2 ).
    setDescription( "Flag that indicates the order of the Ricker to be used o, 1 or 2. Order 2 by default" );

  registerWrapper( viewKeyStruct::outputSeismoTraceString(), &m_outputSeismoTrace ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag that indicates if we write the seismo trace in a file .txt, 0 no output, 1 otherwise" );

  registerWrapper( viewKeyStruct::dtSeismoTraceString(), &m_dtSeismoTrace ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Time step for output pressure at receivers" );

  registerWrapper( viewKeyStruct::indexSeismoTraceString(), &m_indexSeismoTrace ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setDescription( "Count for output pressure at receivers" );

  registerWrapper( viewKeyStruct::forwardString(), &m_forward ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Set to 1 to compute forward propagation" );

  registerWrapper( viewKeyStruct::saveFieldsString(), &m_saveFields ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Set to 1 to save fields during forward and restore them during backward" );

  registerWrapper( viewKeyStruct::shotIndexString(), &m_shotIndex ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Set the current shot for temporary files" );

  registerWrapper( viewKeyStruct::enableLifoString(), &m_enableLifo ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Set to 1 to enable LIFO storage feature" );

  registerWrapper( viewKeyStruct::lifoSizeString(), &m_lifoSize ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( std::numeric_limits< int >::max() ).
    setDescription( "Set the capacity of the lifo storage (should be the total number of buffers to store in the LIFO)" );

  registerWrapper( viewKeyStruct::lifoOnDeviceString(), &m_lifoOnDevice ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -80 ).
    setDescription( "Set the capacity of the lifo device storage (if negative, opposite of percentage of remaining memory)" );

  registerWrapper( viewKeyStruct::lifoOnHostString(), &m_lifoOnHost ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -80 ).
    setDescription( "Set the capacity of the lifo host storage (if negative, opposite of percentage of remaining memory)" );

  registerWrapper( viewKeyStruct::usePMLString(), &m_usePML ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to apply PML" );

  registerWrapper( viewKeyStruct::useDASString(), &m_useDAS ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( WaveSolverUtils::DASType::none ).
    setDescription(
    "Flag to indicate if DAS data will be modeled, and which DAS type to use: \"none\" to deactivate DAS, \"strainIntegration\" for strain integration, \"dipole\" for displacement difference" );

  registerWrapper( viewKeyStruct::linearDASSamplesString(), &m_linearDASSamples ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 5 ).
    setDescription( "Number of sample points to be used for strain integration when integrating the strain for the DAS signal" );

  registerWrapper( viewKeyStruct::linearDASGeometryString(), &m_linearDASGeometry ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Geometry parameters for a linear DAS fiber (dip, azimuth, gauge length)" );

  registerWrapper( viewKeyStruct::linearDASVectorXString(), &m_linearDASVectorX ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "X component of the linear DAS direction vector" );

  registerWrapper( viewKeyStruct::linearDASVectorYString(), &m_linearDASVectorY ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Y component of the linear DAS direction vector" );

  registerWrapper( viewKeyStruct::linearDASVectorZString(), &m_linearDASVectorZ ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Z component of the linear DAS direction vector" );

  registerWrapper( viewKeyStruct::sourceNodeIdsString(), &m_sourceNodeIds ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Indices of the nodes (in the right order) for each source point" );

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstants ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the source for the nodes listed in m_sourceNodeIds" );

  registerWrapper( viewKeyStruct::sourceIsAccessibleString(), &m_sourceIsAccessible ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Flag that indicates whether the source is local to this MPI rank" );

  registerWrapper( viewKeyStruct::receiverNodeIdsString(), &m_receiverNodeIds ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Indices of the nodes (in the right order) for each receiver point" );

  registerWrapper( viewKeyStruct::receiverConstantsString(), &m_receiverConstants ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the receiver for the nodes listed in m_receiverNodeIds" );

  registerWrapper( viewKeyStruct::receiverIsLocalString(), &m_receiverIsLocal ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Flag that indicates whether the receiver is local to this MPI rank" );

  registerWrapper( viewKeyStruct::timestepStabilityLimitString(), &m_timestepStabilityLimit ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Set to 1 to apply a stability limit to the simulation timestep. The timestep used is that given by the CFL condition times the cflFactor parameter." );

  registerWrapper( viewKeyStruct::timeStepString(), &m_timeStep ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setApplyDefaultValue( 1 ).
    setDescription( "TimeStep computed with the power iteration method (if we don't want to compute it, it is initialized with the XML value" );


  registerWrapper( viewKeyStruct::receiverRegionString(), &m_receiverRegion ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Region containing the receivers" );

  registerWrapper( viewKeyStruct::receiverElemString(), &m_receiverElem ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Element containing the receivers" );

  registerWrapper( viewKeyStruct::useTaperString(), &m_useTaper ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to apply taper" );

  registerWrapper( viewKeyStruct::reflectivityCoeffString(), &m_reflectivityCoeff ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.001 ).
    setDescription( "Reflectivity coeff for taper" );

  registerWrapper( viewKeyStruct::xMinTaperString(), &m_xMinTaper ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( {0, 0, 0} ).
    setDescription( "Minimal coordinates for taper (left,bottom,front)" );

  registerWrapper( viewKeyStruct::xMaxTaperString(), &m_xMaxTaper ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( {0, 0, 0} ).
    setDescription( "Maximal coordinates for taper (right,top,back)" );

  registerWrapper( viewKeyStruct::thicknessMinXYZTaperString(), &m_thicknessMinXYZTaper ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( {0, 0, 0} ).
    setDescription( "Size for the taper layer (left,bottom,front)" );

  registerWrapper( viewKeyStruct::thicknessMaxXYZTaperString(), &m_thicknessMaxXYZTaper ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( {0, 0, 0} ).
    setDescription( "Size for the taper layer (left,bottom,front)" );

  registerWrapper( viewKeyStruct::slsReferenceAngularFrequenciesString(), &m_slsReferenceAngularFrequencies ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setApplyDefaultValue( { } ).
    setDescription( "Reference angular frequencies (omega) for the standard-linear-solid (SLS) anelasticity."
                    "The default value is { }, corresponding to no attenuation. An array with the corresponding anelasticity coefficients must be provided." );

  registerWrapper( viewKeyStruct::slsAnelasticityCoefficientsString(), &m_slsAnelasticityCoefficients ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setApplyDefaultValue( { } ).
    setDescription( "Anelasticity coefficients for the standard-linear-solid (SLS) anelasticity."
                    "The default value is { }, corresponding to no attenuation. An array with the corresponding reference frequencies must be provided." );

  registerWrapper( viewKeyStruct::attenuationTypeString(), &m_attenuationType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( WaveSolverUtils::AttenuationType::none ).
    setDescription( "Flag to indicate which attenuation model to use: \"none\" for no attenuation, \"sls\\" " for the standard-linear-solid (SLS) model (Fichtner, 2014)." );

  registerWrapper( viewKeyStruct::sourceWaveletTableNames(), &m_sourceWaveletTableNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the table functions, one for each source, that are used to define the source wavelets. If a list is given, it overrides the Ricker wavelet definitions."
                    "The default value is an empty list, which means that a Ricker wavelet is used everywhere." );

}

WaveSolverBase::~WaveSolverBase()
{
  // TODO Auto-generated destructor stub
}

void WaveSolverBase::reinit()
{
  initializePreSubGroups();
  postInputInitialization();
  initializePostInitialConditionsPreSubGroups();
}

void WaveSolverBase::registerDataOnMesh( Group & meshBodies )
{
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< fields::referencePosition32 >( this->getName() );
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

    nodeManager.getField< fields::referencePosition32 >().resizeDimension< 1 >( X.size( 1 ) );
    arrayView2d< wsCoordType, nodes::REFERENCE_POSITION_USD > const nodeCoords32 = nodeManager.getField< fields::referencePosition32 >();
    for( int i = 0; i < X.size( 0 ); i++ )
    {
      for( int j = 0; j < X.size( 1 ); j++ )
      {
        nodeCoords32[i][j] = X[i][j];
      }
    }


    nodeManager.registerField< fields::taperCoeff >( this->getName());

  } );
}

void WaveSolverBase::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  localIndex const numNodesPerElem = WaveSolverBase::getNumNodesPerElem();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceNodeIds.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceIsAccessible.resize( numSourcesGlobal );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverNodeIds.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverConstants.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverIsLocal.resize( numReceiversGlobal );

  if( m_useSourceWaveletTables )
  {
    FunctionManager const & functionManager = FunctionManager::getInstance();
    m_sourceWaveletTableWrappers.clear();
    for( integer i = 0; i < m_sourceWaveletTableNames.size(); i++ )
    {
      TableFunction const & sourceWaveletTable = functionManager.getGroup< TableFunction >( m_sourceWaveletTableNames[ i ] );
      m_sourceWaveletTableWrappers.emplace_back( sourceWaveletTable.createKernelWrapper() );
    }
  }

}

void WaveSolverBase::postInputInitialization()
{
  SolverBase::postInputInitialization();

  /// set flag PML to one if a PML field is specified in the xml
  /// if counter>1, an error will be thrown as one single PML field is allowed
  integer counter = 0;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  fsManager.forSubGroups< PerfectlyMatchedLayer >( [&] ( PerfectlyMatchedLayer const & )
  {
    counter++;
  } );
  GEOS_THROW_IF( counter > 1,
                 getDataContext() << ": One single PML field specification is allowed",
                 InputError );

  m_usePML = counter;

  if( m_useDAS == WaveSolverUtils::DASType::none && m_linearDASGeometry.size( 0 ) > 0 )
  {
    m_useDAS = WaveSolverUtils::DASType::strainIntegration;
    m_linearDASSamples = 5;
  }

  if( m_useDAS == WaveSolverUtils::DASType::dipole )
  {
    m_linearDASSamples = 2;
  }

  if( m_useDAS != WaveSolverUtils::DASType::none )
  {
    GEOS_LOG_LEVEL_RANK_0( 1, "Modeling linear DAS data is activated" );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "Linear DAS formulation: {}", m_useDAS == WaveSolverUtils::DASType::strainIntegration ? "strain integration" : "displacement difference" ) );

    GEOS_ERROR_IF( m_linearDASGeometry.size( 1 ) != 3,
                   "Invalid number of geometry parameters for the linear DAS fiber. Three parameters are required: dip, azimuth, gauge length" );

    GEOS_ERROR_IF( m_linearDASGeometry.size( 0 ) != m_receiverCoordinates.size( 0 ),
                   "Invalid number of geometry parameters instances for the linear DAS fiber. It should match the number of receivers." );

    m_linearDASVectorX.resize( m_linearDASGeometry.size( 0 ) );
    m_linearDASVectorY.resize( m_linearDASGeometry.size( 0 ) );
    m_linearDASVectorZ.resize( m_linearDASGeometry.size( 0 ) );
    for( int ircv = 0; ircv < m_linearDASGeometry.size( 0 ); ircv++ )
    {
      R1Tensor dasVector = WaveSolverUtils::computeDASVector( m_linearDASGeometry[ ircv ][ 0 ], m_linearDASGeometry[ ircv ][ 1 ] );
      m_linearDASVectorX( ircv ) = dasVector[ 0 ];
      m_linearDASVectorY( ircv ) = dasVector[ 1 ];
      m_linearDASVectorZ( ircv ) = dasVector[ 2 ];
      if( m_useDAS == WaveSolverUtils::DASType::dipole )
      {
        m_linearDASVectorX( ircv ) /= m_linearDASGeometry[ ircv ][ 2 ];
        m_linearDASVectorY( ircv ) /= m_linearDASGeometry[ ircv ][ 2 ];
        m_linearDASVectorZ( ircv ) /= m_linearDASGeometry[ ircv ][ 2 ];
      }
    }
  }

  if( m_attenuationType == WaveSolverUtils::AttenuationType::sls )
  {
    GEOS_THROW_IF( m_slsReferenceAngularFrequencies.size( 0 ) != m_slsAnelasticityCoefficients.size( 0 ),
                   "The number of attenuation anelasticity coefficients for the SLS model must be equal to the number of reference angular frequencies",
                   InputError );
    if( m_slsReferenceAngularFrequencies.size( 0 ) == 0 || m_slsAnelasticityCoefficients.size( 0 ) == 0 )
    {
      m_slsReferenceAngularFrequencies.resize( 1 );
      m_slsReferenceAngularFrequencies[ 0 ] = 2.0 * M_PI * m_timeSourceFrequency;
      m_slsAnelasticityCoefficients.resize( 1 );
      // set the coefficient to a negative value, so that it will be recomputed when the quality factor is available
      m_slsAnelasticityCoefficients[ 0 ] = -1;
    }
  }


  GEOS_THROW_IF( m_sourceCoordinates.size( 0 ) > 0 && m_sourceCoordinates.size( 1 ) != 3,
                 "Invalid number of physical coordinates for the sources",
                 InputError );

  GEOS_THROW_IF( m_receiverCoordinates.size( 0 ) > 0 && m_receiverCoordinates.size( 1 ) != 3,
                 "Invalid number of physical coordinates for the receivers",
                 InputError );

  EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
  real64 const & maxTime = event.getReference< real64 >( EventManager::viewKeyStruct::maxTimeString() );

  if( m_dtSeismoTrace > 0 )
  {
    m_nsamplesSeismoTrace = int( maxTime / m_dtSeismoTrace) + 1;
  }
  else
  {
    m_nsamplesSeismoTrace = 0;
  }

  GEOS_THROW_IF( m_sourceWaveletTableNames.size() > 0 && m_sourceWaveletTableNames.size() != m_sourceCoordinates.size( 0 ),
                 "Invalid number of source wavelet table names. The number of table functions must be equal to the number of sources",
                 InputError );
  m_useSourceWaveletTables = m_sourceWaveletTableNames.size() > 0;

}

real64 WaveSolverBase::solverStep( real64 const & time_n,
                                   real64 const & dt,
                                   integer const cycleNumber,
                                   DomainPartition & domain )
{
  return explicitStep( time_n, dt, cycleNumber, domain );
}

real64 WaveSolverBase::explicitStep( real64 const & time_n,
                                     real64 const & dt,
                                     integer const cycleNumber,
                                     DomainPartition & domain )
{
  if( m_forward )
  {
    return explicitStepForward( time_n, dt, cycleNumber, domain, m_saveFields );
  }
  else
  {
    return explicitStepBackward( time_n, dt, cycleNumber, domain, m_saveFields );
  }
}

localIndex WaveSolverBase::getNumNodesPerElem()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.getGroupPointer< FiniteElementDiscretization >( m_discretizationName );
  GEOS_THROW_IF( feDiscretization == nullptr,
                 getDataContext() << ": FE discretization not found: " << m_discretizationName,
                 InputError );

  localIndex numNodesPerElem = 0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&]( string const &,
                                       MeshLevel const & mesh,
                                       arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    elemManager.forElementRegions( regionNames,
                                   [&] ( localIndex const,
                                         ElementRegionBase const & elemRegion )
    {
      elemRegion.forElementSubRegions( [&]( ElementSubRegionBase const & elementSubRegion )
      {
        finiteElement::FiniteElementBase const &
        fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
        localIndex const numSupportPoints = fe.getNumSupportPoints();
        if( numSupportPoints > numNodesPerElem )
        {
          numNodesPerElem = numSupportPoints;
        }
      } );
    } );


  } );
  return numNodesPerElem;
}

void WaveSolverBase::computeTargetNodeSet( arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
                                           localIndex const subRegionSize,
                                           localIndex const numQuadraturePointsPerElem )
{
  array1d< localIndex > scratch( subRegionSize * numQuadraturePointsPerElem );
  localIndex i = 0;
  for( localIndex e = 0; e < subRegionSize; ++e )
  {
    for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
    {
      scratch[i++] = elemsToNodes( e, q );
    }
  }
  std::ptrdiff_t const numUniqueValues = LvArray::sortedArrayManipulation::makeSortedUnique( scratch.begin(), scratch.end() );

  m_solverTargetNodesSet.insert( scratch.begin(), scratch.begin() + numUniqueValues );
}

void WaveSolverBase::incrementIndexSeismoTrace( real64 const time_n )
{
  while( (m_dtSeismoTrace * m_indexSeismoTrace) <= (time_n + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace )
  {
    m_indexSeismoTrace++;
  }
}

void WaveSolverBase::computeAllSeismoTraces( real64 const time_n,
                                             real64 const dt,
                                             arrayView1d< real32 const > const var_np1,
                                             arrayView1d< real32 const > const var_n,
                                             arrayView2d< real32 > varAtReceivers,
                                             arrayView1d< real32 > coeffs,
                                             bool add )
{
  /*
   * In forward case we compute seismo if time_n + dt is the first time
   * step after the timeSeismo to write.
   *
   *  time_n        timeSeismo    time_n + dt
   *   ---|--------------|-------------|
   *
   * In backward (time_n goes decreasing) case we compute seismo if
   * time_n is the last time step before the timeSeismo to write.
   *
   *  time_n - dt    timeSeismo    time_n
   *   ---|--------------|-------------|
   */

  if( m_nsamplesSeismoTrace == 0 )
    return;
  integer const dir = m_forward ? +1 : -1;
  for( localIndex iSeismo = m_indexSeismoTrace; iSeismo < m_nsamplesSeismoTrace; iSeismo++ )
  {
    real64 const timeSeismo = m_dtSeismoTrace * (m_forward ? iSeismo : (m_nsamplesSeismoTrace - 1) - iSeismo);
    if( dir * timeSeismo > dir * (time_n + epsilonLoc) )
      break;
    WaveSolverUtils::computeSeismoTrace( time_n, dir * dt, timeSeismo, iSeismo, m_receiverNodeIds,
                                         m_receiverConstants, m_receiverIsLocal, var_np1, var_n, varAtReceivers, coeffs, add );
  }
}

void WaveSolverBase::compute2dVariableAllSeismoTraces( localIndex const regionIndex,
                                                       real64 const time_n,
                                                       real64 const dt,
                                                       arrayView2d< real32 const > const var_np1,
                                                       arrayView2d< real32 const > const var_n,
                                                       arrayView2d< real32 > varAtReceivers )
{
  if( m_nsamplesSeismoTrace == 0 )
    return;
  integer const dir = m_forward ? +1 : -1;
  for( localIndex iSeismo = m_indexSeismoTrace; iSeismo < m_nsamplesSeismoTrace; iSeismo++ )
  {
    real64 const timeSeismo = m_dtSeismoTrace * (m_forward ? iSeismo : (m_nsamplesSeismoTrace - 1) - iSeismo);
    if( dir * timeSeismo > dir * (time_n + epsilonLoc))
      break;
    WaveSolverUtils::compute2dVariableSeismoTrace( time_n, dir * dt, regionIndex, m_receiverRegion, timeSeismo, iSeismo, m_receiverElem,
                                                   m_receiverConstants, m_receiverIsLocal, var_np1, var_n, varAtReceivers );
  }
}

bool WaveSolverBase::directoryExists( std::string const & directoryName )
{
  struct stat buffer;
  return stat( directoryName.c_str(), &buffer ) == 0;
}


} /* namespace geos */
