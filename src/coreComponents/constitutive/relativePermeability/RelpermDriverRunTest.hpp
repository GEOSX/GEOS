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

#ifndef GEOSX_RELPERMDRIVERRUNTEST_HPP_
#define GEOSX_RELPERMDRIVERRUNTEST_HPP_

#include "constitutive/relativePermeability/RelpermDriver.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
#include "layouts.hpp"


namespace geosx
{

using namespace constitutive;

//specific to Hysteresis
template< typename RELPERM_TYPE >
std::enable_if_t< std::is_same< TableRelativePermeabilityHysteresis, RELPERM_TYPE >::value, void >
RelpermDriver::runTest( RELPERM_TYPE & relperm,
                        const arrayView2d< real64 > & table )
{
  // get number of phases and components
  integer const numPhases = relperm.numFluidPhases();

  // create kernel wrapper

  typename TableRelativePermeabilityHysteresis::KernelWrapper const kernelWrapper = relperm.createKernelWrapper();

  // set saturation to user specified feed
  // it is more convenient to provide input in molar, so perform molar to mass conversion here

  array2d< real64, compflow::LAYOUT_PHASE > saturationValues;
  if( numPhases > 2 )
  {
    saturationValues.resize( ( m_numSteps + 1 ) * ( m_numSteps + 1 ), numPhases );
  }
  else
  {
    saturationValues.resize( m_numSteps + 1, numPhases );
  }
  using PT = typename RELPERM_TYPE::PhaseType;
  integer const ipWater = relperm.getPhaseOrder()[PT::WATER];
  integer const ipOil = relperm.getPhaseOrder()[PT::OIL];
  integer const ipGas = relperm.getPhaseOrder()[PT::GAS];
  localIndex offset = std::max( std::max( ipOil, ipWater ), std::max( ipOil, ipGas ) ) + 1;

  integer ipWetting = -1, ipNonWetting = -1;
  std::tie( ipWetting, ipNonWetting ) = relperm.phaseIndex();

  for( integer n = 0; n < table.size( 0 ); ++n )
  {
    if( m_numPhases > 2 )
    {
      saturationValues[n][ipWater] = table( n, ipWater + 1 );
      saturationValues[n][ipOil] = table( n, ipOil + 1 );
      saturationValues[n][ipGas] = table( n, ipGas + 1 );
    }
    else//two-phase
    {
      if( ipWater < 0 )
      {
        saturationValues[n][ipOil] = table( n, ipOil + 1 );
        saturationValues[n][ipGas] = table( n, ipGas + 1 );
      }
      else if( ipGas < 0 )
      {
        saturationValues[n][ipWater] = table( n, ipWater + 1 );
        saturationValues[n][ipOil] = table( n, ipOil + 1 );
      }
      else if( ipOil < 0 )
      {
        saturationValues[n][ipWater] = table( n, ipWater + 1 );
        saturationValues[n][ipGas] = table( n, ipGas + 1 );
      }
    }
  }


  arrayView2d< real64 const, compflow::USD_PHASE > const saturation = saturationValues.toViewConst();

  auto const & phaseHasHysteresis = relperm.template getReference< array1d< integer > >( TableRelativePermeabilityHysteresis::viewKeyStruct::phaseHasHysteresisString());

  arrayView2d< real64, compflow::USD_PHASE > phaseMaxHistoricalVolFraction = relperm.template getField< fields::relperm::phaseMaxHistoricalVolFraction >().reference();
  arrayView2d< real64, compflow::USD_PHASE >  phaseMinHistoricalVolFraction = relperm.template getField< fields::relperm::phaseMinHistoricalVolFraction >().reference();

  arrayView1d< real64 > const drainagePhaseMinVolFraction = relperm.template getReference< array1d< real64 > >(
    TableRelativePermeabilityHysteresis::viewKeyStruct::drainagePhaseMinVolumeFractionString());
  arrayView1d< real64 > const drainagePhaseMaxVolFraction = relperm.template getReference< array1d< real64 > >(
    TableRelativePermeabilityHysteresis::viewKeyStruct::drainagePhaseMaxVolumeFractionString());

  //setting for drainage
  {
    if( phaseHasHysteresis[ipNonWetting] )
    {
      phaseMaxHistoricalVolFraction[0][ipNonWetting] = drainagePhaseMaxVolFraction[ipNonWetting];
      std::cout << " new Max NWet Historical " << phaseMaxHistoricalVolFraction[0][ipNonWetting] << std::endl;
    }
    if( phaseHasHysteresis[ipWetting] )
    {
      phaseMinHistoricalVolFraction[0][ipWetting] = drainagePhaseMinVolFraction[ipWetting];
      std::cout << " new Min Wet Historical " << phaseMinHistoricalVolFraction[0][ipWetting] << std::endl;
    }
  }



  forAll< parallelDevicePolicy<> >( saturation.size( 0 ),
                                    [numPhases, kernelWrapper, saturation, table,
                                     offset] GEOSX_HOST_DEVICE ( integer const n )
  {
    // nw phase set max to snw_max to get the imbibition bounding curve
    kernelWrapper.update( 0, 0, saturation[n] );
    for( integer p = 0; p < numPhases; ++p )
    {
      table( n, offset + 1 + p ) = kernelWrapper.relperm()( 0, 0, p );
    }
  } );

  //loop in charge of hysteresis values
  offset += numPhases;

//setting for imbibition
  {
    if( phaseHasHysteresis[ipNonWetting] )
    {
      phaseMaxHistoricalVolFraction[0][ipNonWetting] = drainagePhaseMinVolFraction[ipNonWetting];
      std::cout << " new Max NWet Historical " << phaseMaxHistoricalVolFraction[0][ipNonWetting] << std::endl;
    }
    if( phaseHasHysteresis[ipWetting] )
    {
      phaseMinHistoricalVolFraction[0][ipWetting] = drainagePhaseMaxVolFraction[ipWetting];
      std::cout << " new Min Wet Historical " << phaseMinHistoricalVolFraction[0][ipWetting] << std::endl;
    }
  }



  forAll< parallelDevicePolicy<> >( saturation.size( 0 ),
                                    [numPhases, kernelWrapper, saturation, table,
                                     offset] GEOSX_HOST_DEVICE ( integer const n )
  {
    // nw phase set max to snw_max to get the imbibition bounding curve

    kernelWrapper.update( 0, 0, saturation[n] );
    for( integer p = 0; p < numPhases; ++p )
    {
      table( n, offset + 1 + p ) = kernelWrapper.relperm()( 0, 0, p );
    }
  } );


}

template< typename RELPERM_TYPE >
std::enable_if_t< !std::is_same< TableRelativePermeabilityHysteresis, RELPERM_TYPE >::value, void >
RelpermDriver::runTest( RELPERM_TYPE & relperm,
                        const arrayView2d< real64 > & table )
{
  // get number of phases and components

  integer const numPhases = relperm.numFluidPhases();

  // create kernel wrapper

  typename RELPERM_TYPE::KernelWrapper const kernelWrapper = relperm.createKernelWrapper();

  // set saturation to user specified feed
  // it is more convenient to provide input in molar, so perform molar to mass conversion here

  array2d< real64, compflow::LAYOUT_PHASE > saturationValues;
  if( numPhases > 2 )
  {
    saturationValues.resize(( m_numSteps + 1 ) * ( m_numSteps + 1 ), numPhases );
  }
  else
  {
    saturationValues.resize( m_numSteps + 1, numPhases );
  }
  using PT = typename RELPERM_TYPE::PhaseType;
  integer const ipWater = relperm.getPhaseOrder()[PT::WATER];
  integer const ipOil = relperm.getPhaseOrder()[PT::OIL];
  integer const ipGas = relperm.getPhaseOrder()[PT::GAS];
  const localIndex offset = std::max( std::max( ipOil, ipWater ), std::max( ipOil, ipGas ) ) + 1;

  for( integer n = 0; n < table.size( 0 ); ++n )
  {


    if( m_numPhases > 2 )
    {
      saturationValues[n][ipWater] = table( n, ipWater + 1 );
      saturationValues[n][ipOil] = table( n, ipOil + 1 );
      saturationValues[n][ipGas] = table( n, ipGas + 1 );
    }
    else//two-phase
    {
      if( ipWater < 0 )
      {
        saturationValues[n][ipOil] = table( n, ipOil + 1 );
        saturationValues[n][ipGas] = table( n, ipGas + 1 );
      }
      else if( ipGas < 0 )
      {
        saturationValues[n][ipWater] = table( n, ipWater + 1 );
        saturationValues[n][ipOil] = table( n, ipOil + 1 );
      }
    }

  }

  arrayView2d< real64 const, compflow::USD_PHASE > const saturation = saturationValues.toViewConst();

  // perform relperm update using table (Swet,Snonwet) and save resulting total density, etc.
  // note: column indexing should be kept consistent with output file header below.

  forAll< parallelDevicePolicy<> >( saturation.size( 0 ),
                                    [numPhases, kernelWrapper, saturation, table,
                                     offset] GEOSX_HOST_DEVICE ( integer const n )
  {
    kernelWrapper.update( 0, 0, saturation[n] );
    for( integer p = 0; p < numPhases; ++p )
    {
      table( n, offset + 1 + p ) = kernelWrapper.relperm()( 0, 0, p );
    }
  } );

}


}


#endif //GEOSX_RELPERMDRIVERRUNTEST_HPP_
