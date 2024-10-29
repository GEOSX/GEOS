/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FunctionBase.cpp
 */

#include "FunctionBase.hpp"

namespace geos
{

using namespace dataRepository;



FunctionBase::FunctionBase( const string & name,
                            Group * const parent ):
  Group( name, parent ),
  m_inputVarNames()
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( keys::inputVarNames, &m_inputVarNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Name of fields are input to function." );
}

integer FunctionBase::isFunctionOfTime() const
{
  if( std::find( m_inputVarNames.begin(), m_inputVarNames.end(), "time" ) != m_inputVarNames.end() )
  {
    return 1 + ( m_inputVarNames.size() == 1 );
  }
  return 0;
}

real64_array FunctionBase::evaluateStats( dataRepository::Group const & group,
                                          real64 const time,
                                          SortedArray< localIndex > const & set ) const
{
  localIndex N = set.size();
  real64_array sub( N );
  evaluate( group, time, set.toViewConst(), sub );

  real64_array result( 3 );
  result[0] = 1e10;   // min
  result[1] = 0.0;    // avg
  result[2] = -1e10;  // max
  for( localIndex ii=0; ii<N; ii++ )
  {
    result[0] = std::min( result[0], sub[ii] );
    result[1] += sub[ii];
    result[2] = std::max( result[2], sub[ii] );
  }
  result[1] /= N;

  return result;
}


} // end of namespace geos
