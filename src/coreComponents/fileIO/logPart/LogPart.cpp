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

LogPart::LogPart( string_view logPartTitle )
{
  m_startDesc.m_logPartWidth = std::max((size_t)m_rowMinWidth, logPartTitle.length());
  m_endDesc.m_logPartWidth = std::max((size_t)m_rowMinWidth, logPartTitle.length());

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
  m_rowMinWidth = minWidth;
  m_startDesc.m_logPartWidth = std::max( m_rowMinWidth, m_startDesc.m_logPartWidth );
  m_endDesc.m_logPartWidth = std::max( m_rowMinWidth, m_endDesc.m_logPartWidth );
}

void LogPart::setMaxWidth( size_t const & maxWidth )
{
  m_startDesc.m_logPartMaxWidth = std::min( maxWidth, m_startDesc.m_logPartMaxWidth );
  m_endDesc.m_logPartMaxWidth = std::min( maxWidth, m_endDesc.m_logPartMaxWidth );
}

std::vector< string > splitAndFormatStringByDelimiter( string & description, size_t maxLength, char delimiter = ',' )
{
  std::vector< string > formattedDescription;
  size_t idxBeginCapture = 0;
  size_t idxEndCapture = 0;
  std::vector< std::string > descriptions;
  std::stringstream currLine;
  for( auto strIt = description.begin(); strIt != description.end(); ++strIt )
  {
    if( *strIt == delimiter )
    {
      size_t const dist = std::distance( description.begin(), strIt ) + 1;
      idxEndCapture = dist - idxEndCapture;
      if( dist > maxLength )
      {
        formattedDescription.push_back( currLine.str() );
        currLine.str( "" );
      }
      currLine << description.substr( idxBeginCapture, idxEndCapture );
      idxBeginCapture = dist;
    }
  }
  formattedDescription.push_back( description.substr( idxBeginCapture ) );

  return formattedDescription;
}

void LogPart::formatDescriptions( LogPart::Description & description )
{
  size_t maxNameSize = 1;
  std::vector< string > & descriptionName = description.m_descriptionNames;
  std::vector< std::vector< string > > & m_descriptionsValues = description.m_descriptionsValues;

  for( auto const & name : descriptionName )
  {
    size_t const idx = &name - &(*descriptionName.begin());

    if( !m_descriptionsValues[idx].empty())
    {
      maxNameSize = std::max( maxNameSize, name.size() );
    }
  }

  std::vector< string > & formattedDescriptionLines = description.m_formattedDescriptionLines;
  size_t & logPartWidth = description.m_logPartWidth;
  size_t & logPartMaxWidth = description.m_logPartMaxWidth;
  for( size_t idxName = 0; idxName < descriptionName.size(); idxName++ )
  {
    string & name = descriptionName[idxName];
    std::vector< string > & values = m_descriptionsValues[idxName];

    if( values.empty())
    {
      if( name.size() > logPartMaxWidth )
      {
        std::vector< string > formattedName = splitAndFormatStringByDelimiter( name, logPartMaxWidth );
        for( auto const & format : formattedName )
        {
          formattedDescriptionLines.push_back( format );
        }
      }
      else
      {
        formattedDescriptionLines.push_back( name );
      }
    }
    else
    {
      // +1 for extra space in case of delimiter ":"
      string const spaces = std::string( maxNameSize - name.size() + 1, ' ' );
      string const formattedName = GEOS_FMT( "{}{}", name, spaces );

      string firstValue = values[0];

      if( formattedName.size() + firstValue.size() > logPartMaxWidth )
      {
        std::vector< string > formattedDescription =
          splitAndFormatStringByDelimiter( firstValue, logPartMaxWidth - formattedName.size());
        for( auto const & format : formattedDescription )
        {
          if( &format == &formattedDescription.front())
          {

            formattedDescriptionLines.push_back( GEOS_FMT( "{}{}", formattedName, format ));
          }
          else
          {
            formattedDescriptionLines.push_back(
              GEOS_FMT( "{:>{}}", format, formattedName.size() + format.size() ));
          }
        }
      }
      else
      {
        formattedDescriptionLines.push_back( GEOS_FMT( "{}{}", formattedName, firstValue ));
      }

      logPartWidth = std::max( logPartWidth, formattedDescriptionLines[idxName].size() );

      for( size_t idxValue = 1; idxValue < values.size(); idxValue++ )
      {
        formattedDescriptionLines.push_back(
          GEOS_FMT( "{:>{}}", values[idxValue], formattedName.size() + values[idxValue].size() ));
        logPartWidth = std::max( logPartWidth, formattedDescriptionLines.back().size() );
      }
    }
  }
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

void LogPart::begin( std::ostream & os )
{
  string bottomPart = "";
  if( !m_startDesc.m_descriptionNames.empty())
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
  string const line =  string( m_endDesc.m_logPartWidth, m_borderCharacter );
  if( !m_endDesc.m_descriptionNames.empty() )
  {
    formatDescriptions( m_endDesc );
    topPart = GEOS_FMT( "{}{}\n", buildDescriptionPart( m_endDesc ), line );
  }

  string const bottomPart = GEOS_FMT( "{}{}\n", buildTitlePart( m_endDesc ), line );
  os << GEOS_FMT( "\n{}{}\n", topPart, bottomPart );
}

}
