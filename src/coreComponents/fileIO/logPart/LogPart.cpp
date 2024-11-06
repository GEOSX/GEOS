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
 * @file LogPart.cpp
 */

#include "LogPart.hpp"
#include <algorithm>

namespace geos
{

LogPart::LogPart( string_view logPartTitle ):
  m_logPartTitle( string( logPartTitle ) ),
  m_logPartWidth( std::max((size_t)m_rowMinWidth, logPartTitle.length()) )
{
  m_footerTitle = GEOS_FMT( "End of {}", m_logPartTitle );
}

void LogPart::addDescription( string_view description )
{
  m_beginningDescs.push_back( string( description ) );
}

void LogPart::addEndDescription( string_view description )
{
  m_endDescs.push_back( string( description ) );
}

void LogPart::setMinWidth( integer const & minWidth )
{
  m_rowMinWidth = minWidth;
  m_logPartWidth = std::max( m_rowMinWidth, m_logPartWidth );
}

void LogPart::formatAndInsertDescriptions( std::vector< string > & formattedDescriptions,
                                           string_view name,
                                           std::vector< string > const & descriptionValues )
{
  string const nameFormatted =  GEOS_FMT( "- {}: ", string( name ));

  // format each description with name and the first value
  string const descriptionFormatted = GEOS_FMT( "{}{}", nameFormatted, descriptionValues[0] );

  integer const spacesFromBorder = m_marginBorder * 2 + m_nbBorderChar * 2;
  integer const rowLength = descriptionFormatted.length() + spacesFromBorder;

  formattedDescriptions.push_back( descriptionFormatted );

  m_logPartWidth = std::max( rowLength, m_logPartWidth );

  //format remaining description lines
  for( size_t idxValue = 1; idxValue < descriptionValues.size(); idxValue++ )
  {
    size_t const spaces = descriptionValues[idxValue].length() + nameFormatted.length();
    formattedDescriptions.push_back( GEOS_FMT( "{:>{}}", descriptionValues[idxValue], spaces ) );
  }
}

string LogPart::buildDescriptionPart( std::vector< string > const & formattedDescriptions ) const
{
  std::ostringstream oss;
  for( auto const & description : formattedDescriptions )
  {
    // length of white space to add after the formatted description
    integer const remainingLength = m_logPartWidth - m_nbBorderChar * 2 - m_marginBorder;
    string const borderCharacters = string( m_nbBorderChar, m_borderCharacter );
    oss << borderCharacters;
    oss << GEOS_FMT( "{:<{}}{:<{}}", " ", m_marginBorder, description, remainingLength );
    oss << borderCharacters << '\n';
  }
  return oss.str();
}

string LogPart::buildTitlePart( string_view title ) const
{
  std::ostringstream oss;
  integer const titleRowLength = m_logPartWidth - m_nbBorderChar * 2;
  string const borderCharacters =  string( m_nbBorderChar, m_borderCharacter );
  oss <<  GEOS_FMT( "{}{:^{}}{}\n",
                    borderCharacters,
                    title,
                    titleRowLength,
                    borderCharacters );
  return oss.str();
}

void LogPart::begin( std::ostream & os ) const
{
  string bottomPart;
  if( !m_beginningDescs.empty())
  {
    bottomPart = buildDescriptionPart( m_beginningDescs );
  }

  string const horizontalBorder = string( m_logPartWidth, m_borderCharacter );
  string topPart =  GEOS_FMT( "{}\n{}{}\n", horizontalBorder,
                              buildTitlePart( m_logPartTitle ),
                              horizontalBorder );
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
}

void LogPart::end( std::ostream & os ) const
{
  string topPart;
  string const horizontalBorder =  string( m_logPartWidth, m_borderCharacter );
  if( !m_endDescs.empty() )
  {
    topPart =  GEOS_FMT( "{}{}\n", buildDescriptionPart( m_endDescs ), horizontalBorder );
  }

  string const bottomPart = GEOS_FMT( "{}{}\n", buildTitlePart( m_footerTitle ), horizontalBorder );
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
}

}
