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
 * @file ReactionsBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_REACTIVE_CHEMICALREACTIONS_REACTIONSBASE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_REACTIVE_CHEMICALREACTIONS_REACTIONSBASE_HPP_

#include "dataRepository/ObjectCatalog.hpp"

namespace geos
{

namespace constitutive
{

namespace chemicalReactions
{

class ReactionsBase
{
public:

  ReactionsBase( string const & name, integer const numPrimarySpecies, integer const numSecSpecies );

  virtual ~ReactionsBase() = default;

  constexpr static integer maxNumPrimarySpecies = 12;
  constexpr static integer maxNumSecondarySpecies = 15;

  string const & reactionName() const { return m_name; }

protected:
  /// Name the solubility model
  string m_name;

  integer m_numPrimarySpecies;

  integer m_numSecondarySpecies;

  /// Array storing the name of the components
  string_array m_componentNames;

  /// Array storing the component molar weights
  array1d< real64 > m_componentMolarWeight;

  array1d< real64 > m_log10EqConst;

  array2d< real64 >  m_stoichMatrix;

  array1d< integer > m_chargePrimary;
  array1d< integer > m_chargeSec;

  array1d< real64 >  m_ionSizePrimary;
  array1d< real64 >  m_ionSizeSec;

  real64 m_DebyeHuckelA;
  real64 m_DebyeHuckelB;
  real64 m_WATEQBDot;

  class KernelWrapper
  {
public:

    /**
     * @brief Construct a new Kernel Wrapper object
     *
     * @param log10EqConst
     * @param stoichMatrix
     * @param chargePrimary
     * @param chargeSec
     * @param m_ionSizePrimary
     * @param ionSizeSec
     * @param DebyeHuckelA
     * @param DebyeHuckelB
     * @param WATEQBDot
     */
    KernelWrapper( integer const numPrimarySpecies,
                   integer const numSecondarySpecies,
                   arrayView1d< real64 > const & log10EqConst,
                   arrayView2d< real64 > const & stoichMatrix,
                   arrayView1d< integer > const & chargePrimary,
                   arrayView1d< integer > const & chargeSec,
                   arrayView1d< real64 > const & ionSizePrimary,
                   arrayView1d< real64 > const & ionSizeSec,
                   real64 const DebyeHuckelA,
                   real64 const DebyeHuckelB,
                   real64 const WATEQBDot ):
      m_numPrimarySpecies( numPrimarySpecies ),
      m_numSecondarySpecies( numSecondarySpecies ),
      m_log10EqConst( log10EqConst ),
      m_stoichMatrix( stoichMatrix ),
      m_chargePrimary( chargePrimary ),
      m_chargeSec( chargeSec ),
      m_ionSizePrimary( ionSizePrimary ),
      m_ionSizeSec( ionSizeSec ),
      m_DebyeHuckelA( DebyeHuckelA ),
      m_DebyeHuckelB( DebyeHuckelB ),
      m_WATEQBDot( WATEQBDot )
    {}

protected:

    /**
     * @brief
     *
     * @param temperature
     * @param ionicStrength
     * @param log10PrimaryActCoeff
     * @param dLog10PrimaryActCoeff_dIonicStrength
     * @param log10SecActCoeff
     * @param dLog10SecActCoeff_dIonicStrength
     */
    void computeLog10ActCoefBDotModel( real64 const temperature,
                                       real64 const ionicStrength,
                                       arraySlice1d< real64 > const & log10PrimaryActCoeff,
                                       arraySlice1d< real64 > const & dLog10PrimaryActCoeff_dIonicStrength,
                                       arraySlice1d< real64 > const & log10SecActCoeff,
                                       arraySlice1d< real64 > const & dLog10SecActCoeff_dIonicStrength ) const;
    /**
     * @brief
     *
     * @return
     */
    void computeIonicStrength( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                               arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                               real64 & ionicStrength ) const;


    /// Hard coding the example case - eventually would have to be changed such that it is read from an input file
    integer m_numPrimarySpecies; // Currently not including H2O and O2gas
    integer m_numSecondarySpecies;

    arrayView1d< real64 > m_log10EqConst;
    arrayView2d< real64 > m_stoichMatrix;

    arrayView1d< integer > m_chargePrimary;
    arrayView1d< integer > m_chargeSec;

    arrayView1d< real64 >  m_ionSizePrimary;
    arrayView1d< real64 >  m_ionSizeSec;

    real64 m_DebyeHuckelA;
    real64 m_DebyeHuckelB;
    real64 m_WATEQBDot;
  };

};

} // end namespace chemicalReactions

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_CHEMICALREACTIONS_REACTIONSBASE_HPP_
