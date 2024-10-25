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
 * @file StencilDataCollection.cpp
 */

#include "StencilDataCollection.hpp"

#include "common/Units.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "finiteVolume/TwoPointFluxApproximation.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "common/format/table/TableFormatter.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;


StencilDataCollection::StencilDataCollection( const string & name,
                                              Group * const parent ):
  Base( name, parent )
{
  enableLogLevelInput();
  getWrapperBase( Group::viewKeyStruct::logLevelString() ).
    setDescription( "When higher than 1: Display store events details." );

  registerWrapper( viewKeyStruct::solverNameString(), &m_solverName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the flow solver, to get the permeabilities." );
  registerWrapper( viewKeyStruct::meshNameString(), &m_meshName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of the target " );

  // registerWrapper( viewKeyStruct::connectionDataString(), &m_currentConnData );
  registerWrapper( viewKeyStruct::cellAGlobalIdString(), &m_cellAGlobalId );
  registerWrapper( viewKeyStruct::cellBGlobalIdString(), &m_cellBGlobalId );
  registerWrapper( viewKeyStruct::transmissibilityABString(), &m_transmissibilityAB );
  registerWrapper( viewKeyStruct::transmissibilityBAString(), &m_transmissibilityBA );
}

void StencilDataCollection::postInputInitialization()
{
  Group & problemManager = this->getGroupByPath( "/Problem" );

  { // find targeted solver
    Group & physicsSolverManager = problemManager.getGroup( "Solvers" );

    m_solver = physicsSolverManager.getGroupPointer< FlowSolverBase >( m_solverName );
    GEOS_THROW_IF( m_solver == nullptr,
                   GEOS_FMT( "{}: Could not find flow solver named '{}'.",
                             getDataContext(),
                             m_solverName ),
                   InputError );
  }

  { // find mesh & discretization
//    DomainPartition & domain = problemManager.getDomainPartition();
    DomainPartition & domain = problemManager.getGroup< DomainPartition >( "domain" );

    MeshBody const & meshBody = domain.getMeshBody( m_meshName );
    m_meshLevel = &meshBody.getBaseDiscretization();

    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    try
    {
      string const discretizationName = m_solver->getDiscretizationName();
      m_discretization = &fvManager.getFluxApproximation( discretizationName );
    } catch( BadTypeError const & e )
    {
      // only TPFA is supported for now
      GEOS_ERROR( GEOS_FMT( "{}: target discretization is not supported by {} (for now, only '{}' is).",
                            getDataContext(), catalogName(), TwoPointFluxApproximation::catalogName() ) );
    }
  }
}

void StencilDataCollection::initializePostInitialConditionsPostSubGroups()
{
  // asserting that there is exactly one supported stencil
  integer supportedStencilCount = 0;
  m_discretization->forStencils< CellElementStencilTPFA >( *m_meshLevel, [&]( auto const & stencil )
  {
    globalIndex connCount = stencil.size();
    m_cellAGlobalId.resize( connCount );
    m_cellBGlobalId.resize( connCount );
    m_transmissibilityAB.resize( connCount );
    m_transmissibilityBA.resize( connCount );
    GEOS_LOG_LEVEL_BY_RANK( 1, GEOS_FMT( "{}: initialized {} connection buffer for '{}'.",
                                         getName(), connCount, m_discretization->getName() ) );
    ++supportedStencilCount;
  } );
  GEOS_ERROR_IF( supportedStencilCount == 0, GEOS_FMT( "{}: No compatible discretization was found.", getDataContext() ) );
  GEOS_ERROR_IF( supportedStencilCount > 1, GEOS_FMT( "{}: Multiple discretization was found.", getDataContext() ) );
}


bool StencilDataCollection::execute( real64 const GEOS_UNUSED_PARAM( time_n ),
                                     real64 const GEOS_UNUSED_PARAM( dt ),
                                     integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                     integer const GEOS_UNUSED_PARAM( eventCounter ),
                                     real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                     DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  m_discretization->forStencils< CellElementStencilTPFA >( *m_meshLevel, [&]( auto const & stencil )
  {
    // gather
    auto const stencilWrapper = stencil.createKernelWrapper();
    array1d< KernelConnectionData > const kernelData = gatherConnectionData( stencilWrapper );
    // output
    storeConnectionData( m_discretization->getName(), kernelData.toView() );
  } );

  return false;
}


string showKernelDataExtract( arrayView1d< StencilDataCollection::KernelConnectionData > const & kernelData,
                              globalIndex maxLines = std::numeric_limits< globalIndex >::max() );

class StencilDataCollection::Kernel
{
public:
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using PermeabilityAccessors = StencilMaterialAccessors< PermeabilityBase,
                                                          fields::permeability::permeability >;


  /**
   * @brief launch a kernel to compute the element-element connection data using a given StencilWrapper.
   * @tparam POLICY the Raja host / device launching policy
   * @tparam STENCILWRAPPER_T
   * @param connData The buffer where the connection data will be written
   * @param stencilWrapper the StencilWrapper to use to compute the element-element connection data
   * @param permeability the permeability data that will be used to compute the transmissibility
   */
  template< typename POLICY, typename STENCILWRAPPER_T >
  static void launch( arrayView1d< StencilDataCollection::KernelConnectionData > const & connData,
                      STENCILWRAPPER_T const & stencilWrapper,
                      ElementViewConst< arrayView3d< real64 const > > const & permeability )
  {
    using IndexContainerType = typename STENCILWRAPPER_T::IndexContainerViewConstType;
    IndexContainerType const & elemRegionIndices = stencilWrapper.getElementRegionIndices();
    IndexContainerType const & elemSubRegionIndices = stencilWrapper.getElementSubRegionIndices();
    IndexContainerType const & elementIndices = stencilWrapper.getElementIndices();

    forAll< POLICY >( stencilWrapper.size(), [=] GEOS_HOST_DEVICE ( localIndex const iConn )
    {
      real64 transmissibility[1][2];
      real64 dummy[1][2];

      stencilWrapper.computeWeights( iConn,
                                     permeability,
                                     permeability,      // we don't need derivatives
                                     transmissibility,
                                     dummy );           // we don't need derivatives

      for( localIndex i = 0; i < 2; ++i )
      {
        connData[iConn].m_transmissibility[i] = transmissibility[0][i];
        connData[iConn].m_regionId[i] = elemRegionIndices( iConn, i );
        connData[iConn].m_subRegionId[i] = elemSubRegionIndices( iConn, i );
        connData[iConn].m_elementId[i] = elementIndices( iConn, i );
      }
    } );
    connData.move( LvArray::MemorySpace::host );
  }

private:
  Kernel() = delete;
};

template< typename STENCILWRAPPER_T >
array1d< StencilDataCollection::KernelConnectionData >
StencilDataCollection::gatherConnectionData( STENCILWRAPPER_T const & stencilWrapper ) const
{
  ElementRegionManager const & elemManager = m_meshLevel->getElemManager();

  // allocate a large enough buffer to store all connection data
  array1d< KernelConnectionData > kernelData{ stencilWrapper.size() };
  typename Kernel::PermeabilityAccessors accessor( elemManager, m_solver->getName() );

  Kernel::launch< parallelDevicePolicy<> >( kernelData.toView(), stencilWrapper,
                                            accessor.get< fields::permeability::permeability >() );

  return kernelData;
}


StencilDataCollection::ConnectionData
StencilDataCollection::ConnectionData::fromKernel( KernelConnectionData const & kernelData,
                                                   LocalToGlobalMap const & localToGlobalMap )
{
  ConnectionData conn;
  for( localIndex i = 0; i < 2; ++i )
  {
    conn.m_transmissibility[i] = kernelData.m_transmissibility[i];
    conn.m_globalId[i] = localToGlobalMap[kernelData.m_regionId[i]][kernelData.m_subRegionId[i]][kernelData.m_elementId[i]];
  }
  return conn;
}


string formatKernelDataExtract( arrayView1d< StencilDataCollection::KernelConnectionData > const & kernelData,
                                globalIndex maxLines )
{
  TableData tableData;
  auto kernelIterator = kernelData.begin();
  for( int iConn=0; iConn < maxLines && kernelIterator != kernelData.end(); ++iConn, ++kernelIterator )
  {
    tableData.addRow( GEOS_FMT( "{} / {}", kernelIterator->m_regionId[0], kernelIterator->m_regionId[1] ),
                      GEOS_FMT( "{} / {}", kernelIterator->m_subRegionId[0], kernelIterator->m_subRegionId[1] ),
                      GEOS_FMT( "{} / {}", kernelIterator->m_elementId[0], kernelIterator->m_elementId[1] ),
                      kernelIterator->m_transmissibility[0],
                      kernelIterator->m_transmissibility[1] );
  }
  TableLayout const tableLayout{
    { "regionId A/B", "subRegionId A/B", "elementId A/B", "transmissibilityAB", "transmissibilityBA" },
    GEOS_FMT( "Kernel data (real row count = {})", kernelData.size() )
  };
  TableTextFormatter const tableFormatter{ tableLayout };
  return tableFormatter.toString( tableData );
}

void StencilDataCollection::storeConnectionData( string_view stencilName,
                                                 arrayView1d< KernelConnectionData > const & kernelData )
{
  std::set< ConnectionData > sortedData;

  { // data extraction
    ElementRegionManager const & elemManager = m_meshLevel->getElemManager();
    string const & mapVKStr = ObjectManagerBase::viewKeyStruct::localToGlobalMapString();
    LocalToGlobalMap localToGlobalMap = elemManager.constructArrayViewAccessor< globalIndex, 1 >( mapVKStr );

    std::transform( kernelData.begin(), kernelData.end(),
                    std::inserter( sortedData, sortedData.end() ),
                    [&]( auto const & kernelConn )
    {
      return ConnectionData::fromKernel( kernelConn, localToGlobalMap );
    } );
  }

  { // data storing
    GEOS_ERROR_IF_NE_MSG( size_t( m_cellAGlobalId.size() ), size_t( sortedData.size() ),
                          GEOS_FMT( "{}: Unexpected stencil size!\n{}",
                                    getDataContext(), formatKernelDataExtract( kernelData, 8 ) ) );
    globalIndex i = 0;
    for( ConnectionData const & conn : sortedData )
    {
      m_cellAGlobalId[i] = conn.m_globalId[0];
      m_cellBGlobalId[i] = conn.m_globalId[1];
      m_transmissibilityAB[i] = conn.m_transmissibility[0];
      m_transmissibilityBA[i] = conn.m_transmissibility[1];
      ++i;
    }
  }

  logStoredConnections( stencilName );
}

void StencilDataCollection::logStoredConnections( string_view stencilName )
{
  integer const connCount = MpiWrapper::sum( m_cellAGlobalId.size() );
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}: {} connections stored for '{}'.",
                                      getName(), connCount, stencilName ) );
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        StencilDataCollection,
                        string const &, dataRepository::Group * const )


} /* namespace geos */
