/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "LASWellGenerator.hpp"


namespace geosx
{
using namespace dataRepository;

LASWellGenerator::LASWellGenerator( string const & name, Group * const parent ):
  WellGeneratorBase( name, parent )
{
  registerWrapper(keys::fileName, &m_fileName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setSizedFromParent(0)->
    setDescription("Path to the las file");

  registerWrapper(keys::geometryLogIndexInFile, &m_logIndexToTakeForGeometry, false  )->
    setInputFlag(InputFlags::OPTIONAL)->
    setApplyDefaultValue( -1 )->
    setDescription("Position of the log to take if there are several log sections defined in the LAS file ");
}

LASWellGenerator::~LASWellGenerator()
{
}

void LASWellGenerator::PostProcessInput()
{
}

void LASWellGenerator::GeneratePolyLine()
{
  LASFile lasFile;
  lasFile.Load(  m_fileName );
  if( lasFile.HasLog( "X" ) && lasFile.HasLog( "Y" ) )
  {
    GeneratePolyLineFromXYZ( lasFile );
  }
  else
  {
    GeneratePolyLineFromDepth( lasFile );
  }
  
}

void LASWellGenerator::GeneratePolyLineFromXYZ( LASFile const & lasFile )
{
  auto X = lasFile.GetLog( "X" );
  auto Y = lasFile.GetLog( "Y" );
  auto Z = lasFile.GetLog( "TVDSS" );
  localIndex nbLogEntries = lasFile.LogSize( "X" );
  GEOS_ASSERT( nbLogEntries == lasFile.LogSize( "Y" ) );
  GEOS_ASSERT( nbLogEntries == lasFile.LogSize( "TVDSS" ) );
  m_polyNodeCoords.resize( nbLogEntries );

  auto lasLineXInfo = lasFile.GetLASLines< LASCurveInformationSection >( "X");
  auto lasLineYInfo = lasFile.GetLASLines< LASCurveInformationSection >( "Y" );
  auto lasLineZInfo = lasFile.GetLASLines< LASCurveInformationSection >( "TVDSS");
  GEOS_ERROR_IF(lasLineXInfo.size() > 1, "X curve info is defined in more than one Curve Information Section" );
  GEOS_ERROR_IF(lasLineYInfo.size() > 1, "Y curve info is defined in more than one Curve Information Section" );
  GEOS_ERROR_IF(lasLineZInfo.size() > 1, "TVDSS curve info is defined in more than one Curve Information Section" );
  double factorX = GetFactor( *lasLineXInfo[0] );
  double factorY = GetFactor( *lasLineYInfo[0] );
  double factorZ = GetFactor( *lasLineZInfo[0] );
  m_segmentToPolyNodeMap.resize( nbLogEntries - 1, 2 );
  for( localIndex i = 0; i < nbLogEntries; i++)
  {
    m_polyNodeCoords[i] = { X[i]*factorX, Y[i]*factorY, Z[i]*factorZ };
    if( i < nbLogEntries - 1 )
    {
      m_segmentToPolyNodeMap[i][0] = i;
      m_segmentToPolyNodeMap[i][1] = i +1;
    }
  }
}

void LASWellGenerator::GeneratePolyLineFromDepth( LASFile const & lasFile )
{
  
  auto Xs = lasFile.GetLASLines< LASWellInformationSection >("XCOORD");
  auto Ys = lasFile.GetLASLines< LASWellInformationSection >("YCOORD");
  auto STARTs = lasFile.GetLASLines< LASWellInformationSection >("STRT");
  auto STOPs = lasFile.GetLASLines< LASWellInformationSection >("STOP");
  GEOS_ASSERT( Ys.size() == Xs.size() );
  GEOS_ASSERT( Ys.size() == STARTs.size() );
  GEOS_ASSERT( Ys.size() == STOPs.size() );
  localIndex wellSectionIndex = 0;
  if( Xs.size() > 1 && m_logIndexToTakeForGeometry == -1 )
  {
    GEOS_LOG_RANK_0("Warning : " << this->getName() << " corresponding LAS file has more than 1 log section defined "
                    << "please specify the index of the log section you want to take into account to write the well "
                    << "into the GEOSX data structure. You have to use the keyword " << keys::geometryLogIndexInFile
                    << ". Taking the first one by default.");
  }
  else if( Xs.size() > 1 && m_logIndexToTakeForGeometry != -1 )
  {
    wellSectionIndex = m_logIndexToTakeForGeometry;
  }

  real64 elev = 0.;
  if( lasFile.HasLine<LASWellInformationSection>( "ELEV" ) )
  {
     auto topDepths = lasFile.GetLASLines< LASWellInformationSection >("ELEV");
     elev = topDepths[wellSectionIndex]->GetDataAsReal64();
  }
  real64 topX = Xs[wellSectionIndex]->GetDataAsReal64();
  real64 topY = Ys[wellSectionIndex]->GetDataAsReal64();
  real64 stop = STOPs[wellSectionIndex]->GetDataAsReal64();
  real64 start = STARTs[wellSectionIndex]->GetDataAsReal64();
  R1Tensor top( topX, topY, start - elev );
  R1Tensor bottom( topX, topY, stop - elev );
  if( start < elev ) // Upward positive z axis
  {
    top[2] *= -1;
    bottom[2] *= -1;
  }
  m_polyNodeCoords.resize( 2 );
  m_polyNodeCoords[0] = top;
  m_polyNodeCoords[0] = bottom;
  m_segmentToPolyNodeMap.resize( 1, 2 );
  m_segmentToPolyNodeMap[0][0] = 0;
  m_segmentToPolyNodeMap[0][1] = 1;
  std::cout << top << std::endl;
  std::cout << bottom << std::endl;

}

double LASWellGenerator::GetFactor( LASLine const & lasLine )
{
  double factor = 0.;
  if( stringutilities::ieq( lasLine.GetUnit() , "F" ) ||
      stringutilities::ieq( lasLine.GetUnit() , "ft" ) ||
      stringutilities::ieq( lasLine.GetUnit() , "feets" ) ||
      stringutilities::ieq( lasLine.GetUnit() , "feet" ) )
  {
    factor = 0.3048;
  }
  else if( stringutilities::ieq( lasLine.GetUnit() , "M" ) ||
           stringutilities::ieq( lasLine.GetUnit() , "meters" ) ||
           stringutilities::ieq( lasLine.GetUnit() , "meter" )  ||
           stringutilities::ieq( lasLine.GetUnit() , "metre" ) ||
           stringutilities::ieq( lasLine.GetUnit() , "metres" ) ) 
  {
    factor = 1.;
  }
  else
  {
    GEOS_ERROR( "Unit : " << lasLine.GetUnit() << " is not valid, valid units are FEETS and METERS.");
  }
  return factor;
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, LASWellGenerator, std::string const &, Group * const )

} // namespace
