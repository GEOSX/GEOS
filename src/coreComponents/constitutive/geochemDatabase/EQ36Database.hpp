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
 * @file EQ36Database.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_EQ36DATABASE_HPP
#define GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_EQ36DATABASE_HPP

#include "ThermoDatabaseBase.hpp"

namespace geosx
{

namespace constitutive
{

class EQ36Database : public ThermoDatabaseBase   // placeholder for whatever this class will be derived from
{
public:

  EQ36Database( path const & databaseFileName, string_array const primarySpeciesNames, string_array const secondarySpeciesNames, real64 const temperature );

private:

  // read database and generate species properties, activity coefficients, equilibrium constants and stoichiometric matrix 

  void readThermoDatabaseBase( string_array const & primarySpeciesNames, 
                               string_array const & secondarySpeciesNames, 
                               real64 const temperature ) const;

  void populateActivityCoefficientParameters( real64 const temperature, 
                                              real64 const & DebyeHuckelA, 
                                              real64 const & DebyeHuckelB, 
                                              real64 const & WATEQBDot ) const;

  void populatePrimarySpeciesInformation( string_array const & primarySpeciesNames,
                                          arraySlice1d< real64 > const & MWPrimary,
                                          arraySlice1d< real64 > const & ionSizePrimary,
                                          arraySlice1d< integer > const & chargePrimary ) const;

  void populateSecondarySpeciesInformation( real64 const temperature,
                                            string_array const & secondarySpeciesNames,
                                            arraySlice1d< real64 > const & MWSec,
                                            arraySlice1d< real64 > const & ionSizeSec,
                                            arraySlice1d< integer > const & chargeSec, 
                                            arraySlice2d< real64 > const & stoichMatrix, 
                                            arraySlice1d< real64 > const & log10EqConst ) const;

  void readActivityCoefficientBlock( string const nextBlockString, 
                                     arraySlice1d< real64 > const & readVariable ) const;

  void readAndIgnoreActivityCoefficientBlock( string const nextBlockString ) const;

  void readPrimarySpeciesBlock( string const nextBlockString, 
                                string_array const & primarySpeciesNames, 
                                arraySlice1d< real64 > const & MWPrimary, 
                                arraySlice1d< real64 > const & ionSizePrimary, 
                                arraySlice1d< integer > const & chargePrimary ) const;  

  void readSecondarySpeciesBlock( string const nextBlockString,
                                  unordered_map const primarySpeciesMap,
                                  arraySlice1d< real64 > const temperatureGrid,
                                  unordered_map const & secondarySpeciesMap,
                                  int const & countSecondary, 
                                  arraySlice1d< real64 > const & MWSec,
                                  arraySlice1d< real64 > const & ionSizeSec,
                                  arraySlice1d< integer > const & chargeSec, 
                                  arraySlice2d< real64 > const & stoichMatrix, 
                                  arraySlice1d< real64 > const & log10EqConst ) const;

  void interpolate( arraySlice1d< real64 > const xGrid,
                    arraySlice1d< real64 > const yGrid,
                    real64 const xValue,
                    real64 const & yValue ) const;

  std::ifstream m_EQ36File;
  std::streamsize m_bufferSize = 256;

};

}

}

#endif    // GEOSX_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_EQ36DATABASE_HPP
