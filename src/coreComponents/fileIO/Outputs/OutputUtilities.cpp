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
 * @file OutputUtilities.cpp
 */

#include "OutputUtilities.hpp"

#include "common/logger/Logger.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/NodeManager.hpp"

#include <algorithm>


namespace geos
{

using namespace dataRepository;

class ElementRegionManager;
class NodeManager;

namespace outputUtilities
{

void checkFieldRegistration( ElementRegionManager const & elemManager,
                             NodeManager const & nodeManager,
                             std::set< string > const & fieldNames,
                             string const & outputName )
{
  std::set< string > registeredFieldNames;

  // step 1: loop over all subRegions and scanning cell-based registration
  elemManager.forElementRegions( [&]( ElementRegionBase const & elemRegion )
  {
    elemRegion.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & subRegion )
    {
      for( string const & name : fieldNames )
      {
        if( subRegion.hasWrapper( name ) )
        {
          registeredFieldNames.insert( name );
        }
      }
    } );
  } );

  // step 2: get the fields that have not been found yet on the subRegions
  std::set< string > remainingFieldNames;
  std::set_difference( fieldNames.begin(), fieldNames.end(),
                       registeredFieldNames.begin(), registeredFieldNames.end(),
                       std::inserter( remainingFieldNames, remainingFieldNames.end() ) );

  // step 3: scan node-based registration
  for( string const & name : remainingFieldNames )
  {
    if( nodeManager.hasWrapper( name ) )
    {
      registeredFieldNames.insert( name );
    }
  }

  // step 4: get the fields that are not registered anywhere
  std::set< string > notRegisteredFieldNames;
  std::set_difference( fieldNames.begin(), fieldNames.end(),
                       registeredFieldNames.begin(), registeredFieldNames.end(),
                       std::inserter( notRegisteredFieldNames, notRegisteredFieldNames.end() ) );

  if( !notRegisteredFieldNames.empty() )
  {
    GEOS_LOG_RANK_0( "\n" << outputName << ": Warning! " << notRegisteredFieldNames.size() << " field name specified in `fieldNames` is/are not registered anywhere." );
    for( string const & name : notRegisteredFieldNames )
    {
      GEOS_LOG_RANK_0( outputName << ": Warning! `" << name << "` is not registered anywhere and will not be displayed.\n" );
    }
  }
}

bool isFieldPlotEnabled( PlotLevel const wrapperPlotLevel,
                         PlotLevel const requiredPlotLevel,
                         string const & wrapperName,
                         std::set< string > const & fieldNames,
                         integer const onlyPlotSpecifiedFieldNames )
{
  // check if the plotLevel is sufficiently high for plotting
  bool plotEnabled = ( wrapperPlotLevel <= requiredPlotLevel );

  // if we plot based on plotLevel, then we can return right away if the plotLevel is sufficiently high for plotting
  if( plotEnabled && onlyPlotSpecifiedFieldNames == 0 )
  {
    return plotEnabled;
  }

  // override the logLevel if the fieldNames list was provided
  if( !fieldNames.empty() )
  {
    auto search = fieldNames.find( wrapperName );
    plotEnabled = ( search != fieldNames.end() );
  }
  return plotEnabled;
}


} /* namespace outputUtilities */

} /* namespace geos */
