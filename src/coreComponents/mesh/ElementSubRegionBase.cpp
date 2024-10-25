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
 * @file CellBase.cpp
 */

#include "ElementSubRegionBase.hpp"

namespace geos
{
using namespace dataRepository;

ElementSubRegionBase::ElementSubRegionBase( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_constitutiveModels( groupKeyStruct::constitutiveModelsString(), this ),
  m_numNodesPerElement(),
  m_numEdgesPerElement(),
  m_numFacesPerElement(),
  m_elementCenter(),
  m_elementVolume()
{
  registerGroup( groupKeyStruct::constitutiveModelsString(), &m_constitutiveModels ).
    setSizedFromParent( 1 );

  registerWrapper( viewKeyStruct::numNodesPerElementString(), &m_numNodesPerElement );

  registerWrapper( viewKeyStruct::numEdgesPerElementString(), &m_numEdgesPerElement );

  registerWrapper( viewKeyStruct::numFacesPerElementString(), &m_numFacesPerElement );

  registerWrapper( viewKeyStruct::elementCenterString(), &m_elementCenter ).
    setPlotLevel( PlotLevel::LEVEL_1 ).
    reference().resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::elementVolumeString(), &m_elementVolume ).
    setPlotLevel( PlotLevel::LEVEL_1 );
}

ElementSubRegionBase::~ElementSubRegionBase()
{}

void ElementSubRegionBase::resizePerElementValues( localIndex const newNumNodesPerElement,
                                                   localIndex const newNumEdgesPerElement,
                                                   localIndex const newNumFacesPerElement )
{
  m_numNodesPerElement = newNumNodesPerElement;
  m_numEdgesPerElement = newNumEdgesPerElement;
  m_numFacesPerElement = newNumFacesPerElement;
}


} /* namespace geos */
