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
  m_footerTitle = GEOS_FMT( "{}", m_logPartTitle );
}

void LogPart::addDescription( string const & description )
{
  m_descriptionNames.push_back( description );
  m_descriptionValues.push_back( std::vector< string >() );
}


void LogPart::setMinWidth( size_t const & minWidth )
{
  m_rowMinWidth = minWidth;
  m_logPartWidth = std::max( m_rowMinWidth, m_logPartWidth );
}

void LogPart::formatDescriptions()
{
  size_t maxNameSize = 1;
  for( auto const & name : m_descriptionNames )
  {
    size_t idx = &name - &(*m_descriptionNames.begin());
    if( !m_descriptionValues[idx].empty())
    {
      maxNameSize = std::max( maxNameSize, name.size() );
    }
  }

  size_t const spacesFromBorder = m_borderMargin * 2 + m_nbBorderChar * 2;
  for( size_t idxName = 0; idxName < m_descriptionNames.size(); idxName++ )
  {
    if( m_descriptionValues[idxName].empty())
    {
      m_formattedDescriptions.push_back( m_descriptionNames[idxName] );
    }
    else
    {
      string name = m_descriptionNames[idxName];
      string spaces = string( maxNameSize - name.size(), ' ' );
      string nameFormatted = GEOS_FMT( "{}{} ", name, spaces );

      m_formattedDescriptions.push_back( GEOS_FMT( "{}{}", nameFormatted, m_descriptionValues[idxName][0] ));
      m_logPartWidth = m_formattedDescriptions[idxName].size() + spacesFromBorder;

      for( size_t idxValue = 1; idxValue < m_descriptionValues[idxName].size(); idxValue++ )
      {
        m_formattedDescriptions.push_back( GEOS_FMT( "{:>{}}", m_descriptionValues[idxName][idxValue],
                                                     nameFormatted.size() + m_descriptionValues[idxName][idxValue].size()  ));
        m_logPartWidth = std::max( m_logPartWidth,
                                   m_formattedDescriptions[idxValue].size() + spacesFromBorder );
      }
    }

  }

  m_logPartWidth = std::max( m_rowMinWidth, m_logPartWidth );
}

string LogPart::buildDescriptionPart()
{
  std::ostringstream oss;
  for( auto const & description : m_formattedDescriptions )
  {
    // length of white space to add after the formatted description
    size_t const remainingLength = m_logPartWidth - m_nbBorderChar * 2 - m_borderMargin;
    string const borderCharacters = string( m_nbBorderChar, m_borderCharacter );
    oss << borderCharacters;
    oss << GEOS_FMT( "{:<{}}{:<{}}", " ", m_borderMargin, description, remainingLength );
    oss << borderCharacters << '\n';
  }
  return oss.str();
}

void LogPart::clear()
{
  m_formattedDescriptions.clear();
  m_descriptionNames.clear();
  m_descriptionValues.clear();
}

string LogPart::buildTitlePart()
{
  std::ostringstream oss;
  size_t const titleRowLength = m_logPartWidth - m_nbBorderChar * 2;
  string const borderCharacters =  string( m_nbBorderChar, m_borderCharacter );
  oss <<  GEOS_FMT( "{}{:^{}}{}\n",
                    borderCharacters,
                    m_footerTitle,
                    titleRowLength,
                    borderCharacters );
  return oss.str();
}

void LogPart::begin( std::ostream & os )
{
  string bottomPart = "";
  if( !m_descriptionNames.empty())
  {
    formatDescriptions();
  }
  bottomPart = buildDescriptionPart();

  string const line = string( m_logPartWidth, m_borderCharacter );
  string const topPart =  GEOS_FMT( "{}\n{}{}\n",
                                    line,
                                    buildTitlePart(),
                                    line );
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
  clear();
}

void LogPart::end( std::ostream & os )
{
  string topPart = "";
  string const line =  string( m_logPartWidth, m_borderCharacter );
  if( !m_descriptionNames.empty() )
  {
    formatDescriptions();
    topPart = GEOS_FMT( "{}{}\n", buildDescriptionPart(), line );
  }

  string const bottomPart = GEOS_FMT( "{}{}\n", buildTitlePart(), line );
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
}

}
