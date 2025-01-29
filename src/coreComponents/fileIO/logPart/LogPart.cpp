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
#include "common/format/StringUtilities.hpp"
#include <algorithm>

using namespace geos::stringutilities;
namespace geos
{

LogPart::LogPart( string_view logPartTitle )
{
  m_startDesc.m_title = logPartTitle;
  m_endDesc.m_title = GEOS_FMT( "{}{}", m_prefixEndTitle, logPartTitle );
}

void LogPart::addDescription( string const & description )
{
  m_startDesc.m_descriptionNames.push_back( description );
  m_startDesc.m_descriptionsValues.push_back( std::vector< string >() );
}

void LogPart::addEndDescription( string const & description )
{
  m_endDesc.m_descriptionNames.push_back( description );
  m_endDesc.m_descriptionsValues.push_back( std::vector< string >() );
}


void LogPart::setMinWidth( size_t const & minWidth )
{
  m_startDesc.m_logPartMinWidth = minWidth;
  m_endDesc.m_logPartMinWidth = minWidth;
}

void LogPart::setMaxWidth( size_t const & maxWidth )
{
  m_startDesc.m_logPartMaxWidth = maxWidth;
  m_endDesc.m_logPartMaxWidth = maxWidth;
}

string_view ltrim( string_view s )
{
  std::size_t const first = s.find_first_not_of( " " );
  if( first != string::npos )
  {
    return s.substr( first, ( s.size() - first ) );
  }
  return {};
}

std::vector< string > splitAndFormatStringByDelimiter( string const & description, size_t maxLength )
{
  std::vector< string > formattedDescription;
  size_t startIdx  = 0;
  size_t endIdx = 0;
  size_t captureIdx  = 0;
  size_t spaceIdx = 0;
  while( endIdx < description.size())
  {
    endIdx = description.find( ' ', spaceIdx );
    if( endIdx == std::string::npos )
    {
      if( description.substr( startIdx ).size() >= maxLength )
      {
        formattedDescription.push_back( string( ltrim( description.substr( startIdx, captureIdx ))));
        formattedDescription.push_back( string( ltrim( description.substr( startIdx + captureIdx ))) );
      }
      else
      {
        formattedDescription.push_back( string( ltrim( description.substr( startIdx, endIdx ))));
      }
    }
    else
    {
      size_t partLength = endIdx - startIdx;
      if( partLength >= maxLength )
      {
        formattedDescription.push_back( string( ltrim( description.substr( startIdx, captureIdx ))));
        startIdx = spaceIdx;
        captureIdx = 0;
      }

      spaceIdx = endIdx + 1;
      captureIdx = partLength;
    }
  }
  return formattedDescription;
}

void LogPart::formatDescriptions( LogPart::Description & description )
{
  auto & names = description.m_descriptionNames;
  auto & valuesList = description.m_descriptionsValues;
  size_t & logPartWidth = description.m_logPartWidth;
  size_t & logPartMaxWidth = description.m_logPartMaxWidth;
  size_t & logPartMinWidth = description.m_logPartMinWidth;
  size_t & logPartMaxNameWidth =  description.m_logPartMaxNameWidth;
  std::vector< string > & formattedLines = description.m_formattedDescriptionLines;

  size_t buildingChars = m_nbBorderChar * 2 + m_borderMargin;
  for( size_t idxName = 0; idxName < names.size(); idxName++ )
  {
    string const & name = names[idxName];
    auto const & values = valuesList[idxName];

    if( values.empty())
    {
      if( name.size() > logPartMaxWidth )
      {
        auto formattedName = splitAndFormatStringByDelimiter( name, logPartMaxWidth - buildingChars );
        formattedLines.insert( formattedLines.end(), formattedName.begin(), formattedName.end());
      }
      else
      {
        formattedLines.push_back( name );
      }
    }
    else
    {
      string const spaces = std::string( logPartMaxNameWidth - name.size(), ' ' );
      string const formattedName = GEOS_FMT( "{}{}", name, spaces );
      string const firstValue = values[0];
      size_t const formattedLineWidth = formattedName.size() + firstValue.size() + buildingChars;
      if( formattedLineWidth > logPartMaxWidth )
      {
        auto formattedDescription =
          splitAndFormatStringByDelimiter( firstValue, logPartMaxWidth - formattedName.size() -  buildingChars );
        for( auto const & format : formattedDescription )
        {
          if( &format == &formattedDescription.front())
          {
            formattedLines.push_back( GEOS_FMT( "{}{}", formattedName, format ));
          }
          else
          {
            formattedLines.push_back( GEOS_FMT( "{:>{}}", format, formattedName.size() + format.size() ) );
          }
        }
      }
      else
      {
        formattedLines.push_back( GEOS_FMT( "{}{}", formattedName, firstValue ));
      }

      logPartWidth = std::max( logPartWidth, formattedLines[idxName].size() );

      for( size_t idxValue = 1; idxValue < values.size(); idxValue++ )
      {
        formattedLines.push_back(
          GEOS_FMT( "{:>{}}", values[idxValue], formattedName.size() + values[idxValue].size() ));
        logPartWidth = std::max( logPartWidth, formattedLines.back().size() );
      }
      logPartWidth += m_nbBorderChar * 2;
      if( logPartWidth > logPartMaxWidth )
        logPartWidth = logPartMaxWidth;
    }
  }
  logPartWidth = std::max( logPartWidth, logPartMinWidth );
}

string LogPart::buildDescriptionPart( LogPart::Description const & description )
{
  std::ostringstream oss;
  for( auto const & formattedDescription : description.m_formattedDescriptionLines )
  {
    // length of white space to add after the formatted description
    size_t const remainingLength =  description.m_logPartWidth - m_nbBorderChar * 2 - m_borderMargin;
    string const borderCharacters = string( m_nbBorderChar, m_borderCharacter );
    oss << borderCharacters;
    oss << GEOS_FMT( "{:<{}}{:<{}}", " ", m_borderMargin, formattedDescription, remainingLength );
    oss << borderCharacters << '\n';
  }
  return oss.str();
}

string LogPart::buildTitlePart( LogPart::Description const & description )
{
  std::ostringstream oss;
  size_t const titleRowLength = description.m_logPartWidth - m_nbBorderChar * 2;
  string const borderCharacters =  string( m_nbBorderChar, m_borderCharacter );
  oss << GEOS_FMT( "{}{:^{}}{}\n",
                   borderCharacters,
                   description.m_title,
                   titleRowLength,
                   borderCharacters );
  return oss.str();
}

void LogPart::computeInitialLogWidth( LogPart::Description & description,
                                      std::vector< string > const & descriptionNames,
                                      std::vector< std::vector< string > > m_descriptionsValues )
{

  if( !descriptionNames.empty() )
  {
    size_t maxStringSize = 0;
    for( size_t idxDescription = 0; idxDescription  < descriptionNames.size(); idxDescription++ )
    {
      if( !m_descriptionsValues[idxDescription].empty())
      {
        maxStringSize = std::max( maxStringSize, descriptionNames[idxDescription].size() );
      }
    }

    description.m_logPartMaxNameWidth = maxStringSize;
  }
}

void LogPart::begin( std::ostream & os )
{
  string bottomPart = "";
  auto & descriptionNames = m_startDesc.m_descriptionNames;

  computeInitialLogWidth( m_startDesc, descriptionNames, m_startDesc.m_descriptionsValues );

  if( !descriptionNames.empty())
  {
    formatDescriptions( m_startDesc );
  }
  bottomPart = buildDescriptionPart( m_startDesc );

  string const line = string( m_startDesc.m_logPartWidth, m_borderCharacter );
  string const topPart =  GEOS_FMT( "{}\n{}{}\n",
                                    line,
                                    buildTitlePart( m_startDesc ),
                                    line );
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
}

void LogPart::end( std::ostream & os )
{
  string topPart = "";
  auto & descriptionNames = m_endDesc.m_descriptionNames;

  computeInitialLogWidth( m_endDesc, descriptionNames, m_endDesc.m_descriptionsValues );

  formatDescriptions( m_endDesc );
  string const line =  string( m_endDesc.m_logPartWidth, m_borderCharacter );
  if( !descriptionNames.empty() )
  {
    topPart = GEOS_FMT( "{}{}\n", buildDescriptionPart( m_endDesc ), line );
  }

  string const bottomPart = GEOS_FMT( "{}{}\n", buildTitlePart( m_endDesc ), line );
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
}

}
