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
 * @file EQ36Database.cpp
 */

#include "EQ36Database.hpp"

using namespace std;

namespace geosx
{
using namespace stringutilities;

namespace constitutive
{
EQ36Database::EQ36Database( path const & databaseFileName, string_array const primarySpeciesNames, string_array const secondarySpeciesNames, real64 const temperature ):
  ThermoDatabaseBase( databaseFileName )
{
  readThermoDatabaseBase( primarySpeciesNames, secondarySpeciesNames, temperature );
}

void EQ36Database::readThermoDatabaseBase( string_array const & primarySpeciesNames, 
                                           string_array const & secondarySpeciesNames,
                                           real64 const temperature )
{
  // Read the database file
  std::ifstream m_EQ36File ( m_databaseFileName );

  // Populate activity coefficient related properties
  real64 DebyeHuckelA;
  real64 DebyeHuckelB;
  real64 WATEQBDot;
  populateActivityCoefficientParameters( temperature, DebyeHuckelA, DebyeHuckelB, WATEQBDot );


  // Populate primary species related properties
  integer numPrimarySpecies = primarySpeciesNames.size();
  array1d< real64 > MWPrimary;
  array1d< real64 > ionSizePrimary;
  array1d< int > chargePrimary;
  MWPrimary.resize( numPrimarySpecies );
  ionSizePrimary.resize( numPrimarySpecies );
  chargePrimary.resize( numPrimarySpecies );
  populatePrimarySpeciesParameters( primarySpeciesNames, MWPrimary, ionSizePrimary, chargePrimary  );


  // Populate secondary species related properties

  // flag to determine is user has provided a list of secondary species to track
  // Code is currently not set up to work when user does not list secondary species due to difficulty in working with unknown sizes
  // To allow for automatic population of secondary species based on the database a function similar to readSeconarySpeciesBlock would 
  // have to be written that just counts the number of secondary species and creates the secondarySpeciesNames variable
  bool inputSecondarySpeciesFlag = secondarySpeciesNames.size() == 0 ? 1 : 0; 
  GEOSX_ERROR_IF( inputSecondarySpeciesFlag, "Currently user is expected to supply secondary species names as well." );

  integer numSecondarySpecies = secondarySpeciesNames.size();
  array1d< real64 > MWSec;
  array1d< real64 > ionSizeSec;
  array1d< int > chargeSec;
  array2d< real64 > stoichMatrix;
  array1d< real64 > log10EqConst;
  MWSec.resize( numSecondarySpecies );
  ionSizeSec.resize( numSecondarySpecies );
  chargeSec.resize( numSecondarySpecies );
  stoichMatrix.resize( numSecondarySpecies, numPrimarySpecies );
  log10EqConst.resize( numSecondarySpecies );

  populateSecondarySpeciesInformation( temperature, secondarySpeciesNames, MWSec, ionSizeSec, chargeSec, stoichMatrix, log10EqConst );

  is.close();
}


void EQ36Database::populateActivityCoefficientParameters( real64 const temperature, 
                                                          real64 const & DebyeHuckelA, 
                                                          real64 const & DebyeHuckelB, 
                                                          real64 const & WATEQBDot ) const
{

  // This function reads the activity coefficient related blocks and populates the DebyeHuckel and B-Dot variables. 

  // read and ignore till you reach the temeprature grid block
  nextBlockString = "Temperature grid";
  readAndIgnoreActivityCoefficientBlock( nextBlockString );

  // Ask Matteo how to size the variables used below. 
  // Is it ok to declare it as a bigger array and then resize it to a smaller one?

  // read and store temperature information
  nextBlockString = "Pressure grid";
  readActivityCoefficientBlock( nextBlockString, temperatureGrid );
 
  // read and ignore pressure and pressure envelope information
  // As of now, we don't make use of pressure and hence I am not storing it
  nextBlockString = "Pressure envelope";
  readAndIgnoreActivityCoefficientBlock( nextBlockString );
  nextBlockString = "Debye-Huckel A_gamma";
  readAndIgnoreActivityCoefficientBlock( nextBlockString );
 
  // read and store Debye-Huckel A information
  nextBlockString = "Debye-Huckel A_H";
  readActivityCoefficientBlock( nextBlockString, DebyeHuckelAGrid );
 
  // read and ignore Debye-Huckel A_H information
  nextBlockString = "Debye-Huckel B_gamma";
  readAndIgnoreActivityCoefficientBlock( nextBlockString );
 
  // read and store Debye-Huckel B information and convert it in the correct units
  nextBlockString = "Debye-Huckel B_H";
  readActivityCoefficientBlock( nextBlockString, DebyeHuckelBGrid );
  DebyeHuckelBGrid = DebyeHuckelBGrid*1e8;
 
  // read and ignore Debye-Huckel B_H information
  nextBlockString = "B-dot";
  readAndIgnoreActivityCoefficientBlock( nextBlockString );
 
  // read and store WATEQB-Dot information
  nextBlockString = "B-dot_H";
  readActivityCoefficientBlock( nextBlockString, WATEQBDotGrid );

  // if the size of the variables to be interpolated are not the same, print an error statement
  GEOSX_ERROR_IF( temperatureGrid.size() != DebyeHuckelAGrid.size() || temperatureGrid.size() != DebyeHuckelBGrid.size() || 
                  temperatureGrid.size() != WATEQBDotGrid.size(), "The interpolation table for temperature grid and activity 
                  coefficient parameters have different sizes" );

  // linear interpolation for the parameter values
  interpolate(temperatureGrid, DebyeHuckelAGrid, temperature, DebyeHuckelA)
  interpolate(temperatureGrid, DebyeHuckelBGrid, temperature, DebyeHuckelB)
  interpolate(temperatureGrid, WATEQBDotGrid, temperature, WATEQBDot)

}

void EQ36Database::populatePrimarySpeciesInformation( string_array const & primarySpeciesNames,
                                                        arraySlice1d< real64 > const & MWPrimary,
                                                        arraySlice1d< real64 > const & ionSizePrimary,
                                                        arraySlice1d< integer > const & chargePrimary ) const
{

  // This function reads the primary species related blocks and populates the molecular weight, ion size and charge variables. 

  // read and ignore everything till you reach the basis species block
  nextBlockString = "basis species";
  readAndIgnoreActivityCoefficientBlock( nextBlockString ); 

  // read and store primary species information
  nextBlockString = "auxiliary basis species";
  readPrimarySpeciesBlock( nextBlockString, primarySpeciesNames, MWPrimary, ionSizePrimary, chargePrimary );

}


void EQ36Database::populateSecondarySpeciesInformation( real64 const temperature,
                                                        string_array const & secondarySpeciesNames,
                                                        arraySlice1d< real64 > const & MWSec,
                                                        arraySlice1d< real64 > const & ionSizeSec,
                                                        arraySlice1d< integer > const & chargeSec, 
                                                        arraySlice2d< real64 > const & stoichMatrix, 
                                                        arraySlice1d< real64 > const & log10EqConst ) const


{

  // This function reads the secondary species related blocks and populates the molculear weight, ion size, charge, stoichiometric
  // matrix and equilibrium constant variables.  

  // create a map to get the index numbers of the secondary species
  unordered_map< string, int > secondarySpeciesMap;
  for( int ic = 0; ic < secondarySpeciesNames.size(); ic++ )
    secondarySpeciesMap[ secondarySpeciesNames[ic] ] = ic;

  int countSecondary=0

  // read auxiliary basis species information
  nextBlockString = "aqueous species";
  readSecondarySpeciesBlock( nextBlockString, primarySpeciesMap, temperatureGrid, secondarySpeciesMap, countSecondary, MWSec, ionSizeSec, chargeSec, stoichMatrix, log10EqConst );

  // read auxiliary aqueous species information
  nextBlockString = "solids";
  readSecondarySpeciesBlock( nextBlockString, primarySpeciesMap, temperatureGrid, secondarySpeciesMap, countSecondary, MWSec, ionSizeSec, chargeSec, stoichMatrix, log10EqConst );

  // read solid species information
  nextBlockString = "liquids";
  readSecondarySpeciesBlock( nextBlockString, primarySpeciesMap, temperatureGrid, secondarySpeciesMap, countSecondary, MWSec, ionSizeSec, chargeSec, stoichMatrix, log10EqConst );

  // read liquid species information
  nextBlockString = "gases";
  readSecondarySpeciesBlock( nextBlockString, primarySpeciesMap, temperatureGrid, secondarySpeciesMap, countSecondary, MWSec, ionSizeSec, chargeSec, stoichMatrix, log10EqConst );

  // read gases species information
  nextBlockString = "solid solutions";
  readSecondarySpeciesBlock( nextBlockString, primarySpeciesMap, temperatureGrid, secondarySpeciesMap, countSecondary, MWSec, ionSizeSec, chargeSec, stoichMatrix, log10EqConst );

  GEOSX_ERROR_IF( countSecondary != secondarySpeciesNames.size(), "Could not find all secondary species." );

}

void EQ36Database::readActivityCoefficientBlock( string const nextBlockString,
                                                 arraySlice1d< real64 > const & readVariable ) const

{

  // This function reads a block from the EQ3/6 Database file related to the B-Dot activity coefficient model and returns the read entries.

  char fileLine[m_bufferSize];
  while( m_EQ36File.getline( fileLine, m_bufferSize ))         //read the line
  {
    std::string fileLineString ( fileLine );
    auto found = fileLineString.find( nextBlockString );   // can you find the next block string string in the line
    if( found!=std::string::npos )                         // if yes, you have reached the next block and should exit
      break;

    string_array lineEntries = Tokenize( fileLineString, " " );      // separate the different numeric enteries and store them
    for( int i = 0; i < lineEntries.size(); i++ )
    {
      // Ask Matteo if this is acceptable as readVariable hasn't been resized anywhere
      readVariable.emplace_back( std::stod( lineEntries[i] ));       
    }
  }

}

void EQ36Database::readAndIgnoreActivityCoefficientBlock( string const nextBlockString ) const

{

  // This function reads a block from the EQ3/6 Database file and ignores it as it isn't currently useful

  char fileLine[m_bufferSize];
  while( m_EQ36File.getline( fileLine, m_bufferSize ))         //read the line
  {
    std::string fileLineString ( fileLine );
    auto found = fileLineString.find( nextBlockString );   // can you find the next block string string in the line
    if( found!=std::string::npos )                         // if yes, you have reached the next block and should exit
      break;

  }

}

void EQ36Database::readPrimarySpeciesBlock( string const nextBlockString,
                                            string_array const & primarySpeciesNames,
                                            arraySlice1d< real64 > const & MWPrimary,
                                            arraySlice1d< real64 > const & ionSizePrimary,
                                            arraySlice1d< integer > const & chargePrimary ) const
{

  // This function reads the primary species block and stores the read entries

  // create a map to get the index numbers of the primary species
  unordered_map< string, int > primarySpeciesMap;
  for( int ic = 0; ic < primarySpeciesNames.size(); ic++ )
  {
    primarySpeciesMap[primarySpeciesNames[ic] ] = ic;
  }

  // define variables to store things
  real64 MW;
  real64 ionSize;
  integer charge;
  string speciesName = " ";

  int countPrimary = 0;
  while( m_EQ36File.getline( fileLine, m_bufferSize ))          // read the line
  {
    std::string fileLineString( fileLine );
    {
      auto found = fileLineString.find( "+-------" );       // if one finds "+-----" in the line, it means that the next line is either the start of a new block or new primary species
      if( found!=std::string::npos )                        // Either you have all the information for the previous species or you are just starting
      {
        auto it = primarySpeciesMap.find( speciesName );          
        if( it != primarySpeciesMap.end() )
        {
          // Ignoring speciesType information since it doesn't make a difference at this point
          int speciesIndex = primarySpeciesMap[ speciesName ]; 
          MWPrimary[speciesIndex] = MW;
          ionSizePrimary[speciesIndex] = ionSize;
          chargePrimary[speciesIndex] = charge;
          countPrimary++;
        }

        m_EQ36File.getline( fileLine, m_bufferSize );          // read next line
        speciesName = fileLine;

        auto found2 = speciesName.find( nextBlockString );  // if the next block has started, exit the loop
        if( found2 != std::string::npos )
        {
          break;
        }
      }
    }

    // molecular weight section
    {
      auto found = fileLineString.find( "mol.wt." );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        MW = std::stod( lineEntries[3] ) * 0.001;
      }
    }

    // ion size section
    {
      auto found = fileLineString.find( "DHazero" );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        ionSize = std::stod( lineEntries[3] ) * 1e-8;
      }
    }
    // charge section
    {
      auto found = fileLineString.find( "charge" );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        charge = std::stoi( lineEntries[2] );
      }
    }
  }

  GEOSX_ERROR_IF( countPrimary != primarySpeciesNames.size(), "Could not find all input primary species." );

}

void EQ36Database::readSecondarySpeciesBlock( string const nextBlockString,
                                              unordered_map const primarySpeciesMap,
                                              arraySlice1d< real64 > const temperatureGrid,
                                              unordered_map const & secondarySpeciesMap,
                                              int const & countSecondary, 
                                              arraySlice1d< real64 > const & MWSec,
                                              arraySlice1d< real64 > const & ionSizeSec,
                                              arraySlice1d< integer > const & chargeSec, 
                                              arraySlice2d< real64 > const & stoichMatrix, 
                                              arraySlice1d< real64 > const & log10EqConst ) const
{

  // This function reads a block from the EQ3/6 Database file related to the secondary species and returns the read entries.
  // This function has a high chance of not functioning as intended (due to higher complexity) and needs to be tested for all types of scenarios. 

  // variables to store read information
  string_array speciesNames;
  array1d< integer > speciesIndices;
  array1d< real64 > stoichs;
  array1d< real64 > logKs;

  bool primarySpeciesFound = 0;

  while( m_EQ36File.getline( fileLine, m_bufferSize ))          // read from file
  {
    std::string fileLineString( fileLine );
    {
      auto found = fileLineString.find( "+-------" );       // If +----- is found, it means you are at the beginning of the next species or the beginning of a new block
      if( found!=std::string::npos )
      {
        // record these entries only when all the products in the reaction associated with this species are listed as a primary species in the input file
        // Have to update how that is determined in the bottom. Will do in next round
        if( primarySpeciesFound )
        {
          // Find if the species is listed as a secondary species 
          auto iter = secondarySpeciesMap.find( speciesName );
          if( iter != secondarySpeciesMap.end() )
          {

            int secSpeciesIndex = secondarySpeciesMap[speciesName]

            // Ignoring speciesType information since I don't know what difference does it make
            MWSec[secSpeciesIndex] = MW;
            ionSizeSec[secSpeciesIndex] = ionSize;
            chargeSec[secSpeciesIndex] = charge;
            stoichMatrix[secSpeciesIndex][:] = stoichs;                   // rectify syntax

            GEOSX_ERROR_IF( logKs.size() != temperatureGrid.size(), "The interpolation table for temperature grid and equilibrium constant 
                            have different sizes" );
            real64 tempVariable;
            interpolate( temperatureGrid, logKs, temperature, tempVariable )
            log10EqConst[secSpeciesIndex] = tempVariable

            countSecondary++;
          }

          stoichs.clear();
          logKs.clear();
        }

        m_EQ36File.getline( fileLine, m_bufferSize );                       // read the next line
        std::string fileLineString2( fileLine );
        string_array lineEntries2 = Tokenize( fileLineString2, " " );
        speciesName = lineEntries2[0];

        auto found2 = lineEntries2.find( nextBlockString );           // if reached next block, exit the loop
        if( found2 != std::string::npos )
        {
          break;
        }
      }
    }

    // Molecular weight section
    {
      auto found = fileLineString.find( "mol.wt." );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        MW = std::stod( lineEntries[3] ) * 0.001;
      }
    }

    // Ion size section
    {
      auto found = fileLineString.find( "DHazero" );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        ionSize = std::stod( lineEntries[3] ) * 1e-8;
      }
    }

    // charge section
    {
      auto found = fileLineString.find( "charge" );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        charge = std::stoi( lineEntries[2] );
      }
    }
    // stoichiometric section
    {
      auto found = fileLineString.find( "aqueous dissociation reaction" );
      if( found != std::string::npos )
      {
        string_array lineEntries = Tokenize( fileLineString, " " );
        integer numRxnSpecies =  std::stoi( lineEntries[0] );       // number of species that participate in the reaction
        reactionPrimarySpeciesNames.clear();    
        reactionPrimarySpeciesIndices.clear();    
        stoichs.clear();
        reactionPrimarySpeciesNames.resize( numRxnSpecies-1 );
        reactionPrimarySpeciesIndices.resize( numRxnSpecies-1 );
        stoichs.resize( numRxnSpecies-1 );

        int reactionSpeciesCount = 0;
        while( m_EQ36File.getline( fileLine, m_bufferSize ))
        {
          std::string fileLineString2( fileLine );
          auto found2 = fileLineString2.find( "* Log K" );      // reached the end of this block so exit loop
          if( found2 != std::string::npos )
            break;

          string_array lineEntries2 = Tokenize( fileLineString2, " " );
          for( int i = 0; i < lineEntries2.size(); i = i+2 )
          {
            GEOSX_ERROR_IF( reactionSpeciesCount == 0 & lineEntries2[1] != speciesName, "First species in the reaction is expected to be the secondary species" );
            GEOSX_ERROR_IF( reactionSpeciesCount == 0 & lineEntries2[0] != -1, "First species in the reaction is expected to have a stoichiometry of -1" );
            reactionSpeciesCount++;
            if( reactionSpeciesCount > 1 )
            {
              auto it = primarySpeciesMap.find( lineEntries2[i+1] );
              // Note that no special treatment for H2O or O2 gas has been implemented. 
              // If these are present in one's problem, they need to be listed as input primary species. 
              if ( it != primarySpeciesMap.end() )
              {
                reactionPrimarySpeciesNames[reactionSpeciesCount-2] = lineEntries2[i+1];
                stoichs[reactionSpeciesCount-2] = lineEntries2[i];
                reactionPrimarySpeciesIndices[reactionSpeciesCount-2] = primarySpeciesMap[ lineEntries2[i+1] ];
                primarySpeciesFound = 1;
              }
              else
              {
                primarySpeciesFound = 0;
                break
              }
            }
          }
        }

        //GEOSX_ERROR_IF( reactionSpeciesCount != numRxnSpecies, "Internal error when reading database" );    // This needs to go somewhere else

        // If even one of the species in the reaction is not listed as a primary species in the input file, do not include this secondary species
        // Else, read the equilibrium constant section.         
        if( primarySpeciesFound )
        {
          logKs.clear();
          while( m_EQ36File.getline( fileLine, m_bufferSize ))
          {
            std::string fileLineString2( fileLine );
            auto found2 = fileLineString2.find( "*" );          // if reached next section, exit the loop
            if( found2 != std::string::npos )
              break;

            string_array lineEntries2 = Tokenize( fileLineString2, " " );
            for( int i = 0; i < lineEntries2.size(); i++ )
              logKs.emplace_back( std::stod( lineEntries2[i] ));      // Check with Matteo what to do about this since we don't know the size of the vector
          }
        }
      }
    }
  }
}


void EQ36Database::interpolate( arraySlice1d< real64 > const xGrid,
                                arraySlice1d< real64 > const yGrid,
                                real64 const xValue,
                                real64 const & yValue ) const
{
  // Find the nearest x Value
  int loopIndex = 0;
  while ( xValue < xGrid[loopIndex] & loopIndex < xGrid.size() )
    loopIndex++;
  GEOSX_ERROR_IF( loopIndex == 0 || loopIndex == xGrid.size(), "Input temperature not in permitted range." );   // So far x is always temperature

  // Interpolate 
  yValue = yGrid[loopIndex-1] + (yGrid[loopIndex] - yGrid[loopIndex-1]) / (xGrid[loopIndex] - xGrid[loopIndex-1]) * (xValue - xGrid[loopIndex-1])

}


// Not sure what is this for
REGISTER_CATALOG_ENTRY( ThermoDatabaseBase,
                        EQ36Database,
                        const Path &,
                        const string_array &,
                        const string_array & )

}

}
