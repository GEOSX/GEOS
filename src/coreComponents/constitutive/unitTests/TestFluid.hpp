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
 * @file TestFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUID_HPP_
#define GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUID_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"

namespace geos
{

namespace testing
{

struct Fluid
{
  static constexpr integer NC =  11;

  static constexpr integer H2O =  0;
  static constexpr integer CO2 =  1;
  static constexpr integer N2  =  2;
  static constexpr integer H2S =  3;
  static constexpr integer C1  =  4;
  static constexpr integer C2  =  5;
  static constexpr integer C3  =  6;
  static constexpr integer C4  =  7;
  static constexpr integer C5  =  8;
  static constexpr integer C8  =  9;
  static constexpr integer C10 = 10;

  static constexpr integer Pc = 0;    // Critical pressure
  static constexpr integer Tc = 1;    // Critical temperature
  static constexpr integer Vc = 2;    // Critical colume
  static constexpr integer Ac = 3;    // Accentric factor
  static constexpr integer Mw = 4;    // Molecular weight
  static constexpr integer Vs = 5;    // Volume shift

  static std::array< real64, 66 > data;
};

template< int NC >
using Feed = std::array< real64, NC >;

template< int NC >
class TestFluid
{
public:
  ~TestFluid() = default;

  static std::unique_ptr< TestFluid< NC > > create( std::array< integer, NC > const & components )
  {
    std::unique_ptr< TestFluid< NC > > testFluid( new TestFluid() );
    const std::unordered_map< integer, string > componentNames = {
      {Fluid::H2O, "H2O"}, {Fluid::CO2, "CO2"}, {Fluid::N2, "N2"}, {Fluid::H2S, "H2S"},
      {Fluid::C1, "CH4"}, {Fluid::C2, "C2H6"}, {Fluid::C3, "C3H8"}, {Fluid::C4, "C4H10"},
      {Fluid::C5, "C5H12"}, {Fluid::C8, "C8H18"}, {Fluid::C10, "C10+"},
    };
    for( integer const ic : components )
    {
      testFluid->componentNames.emplace_back( componentNames.at( ic ) );
    }
    createArray( testFluid->criticalPressure, components, Fluid::Pc, Fluid::data );
    createArray( testFluid->criticalTemperature, components, Fluid::Tc, Fluid::data );
    createArray( testFluid->criticalVolume, components, Fluid::Vc, Fluid::data );
    createArray( testFluid->acentricFactor, components, Fluid::Ac, Fluid::data );
    createArray( testFluid->molecularWeight, components, Fluid::Mw, Fluid::data );
    createArray( testFluid->volumeShift, components, Fluid::Vs, Fluid::data );
    testFluid->binaryCoeff.resize( NC, NC );
    return testFluid;
  }

  template< typename LIST >
  void setBinaryCoefficients( LIST const & bics )
  {
    auto bic = bics.begin();
    for( int i = 0; i < NC; ++i )
    {
      for( int j = 0; j < i; ++j )
      {
        binaryCoeff( i, j ) = *bic++;
        binaryCoeff( j, i ) = binaryCoeff( i, j );
      }
      binaryCoeff( i, i ) = 0.0;
    }
  }

  constitutive::compositional::ComponentProperties & getComponentProperties()
  {
    if( !m_component_properties )
    {
      m_component_properties = std::make_unique< constitutive::compositional::ComponentProperties >(
        componentNames,
        molecularWeight );
      createArray( m_component_properties->m_componentCriticalPressure, criticalPressure );
      createArray( m_component_properties->m_componentCriticalTemperature, criticalTemperature );
      createArray( m_component_properties->m_componentAcentricFactor, acentricFactor );
      createArray( m_component_properties->m_componentVolumeShift, volumeShift );
      m_component_properties->m_componentBinaryCoeff.resize( NC, NC );
      for( integer ic = 0; ic < NC; ++ic )
      {
        for( integer jc = 0; jc < NC; ++jc )
        {
          m_component_properties->m_componentBinaryCoeff( ic, jc ) = binaryCoeff( ic, jc );
        }
      }
    }
    return *m_component_properties;
  }

  constitutive::compositional::ComponentProperties::KernelWrapper createKernelWrapper()
  {
    return getComponentProperties().createKernelWrapper();
  }

private:
  TestFluid() = default;

public:
  string_array componentNames;
  array1d< real64 > criticalPressure;
  array1d< real64 > criticalTemperature;
  array1d< real64 > criticalVolume;
  array1d< real64 > acentricFactor;
  array1d< real64 > molecularWeight;
  array1d< real64 > volumeShift;
  array2d< real64 > binaryCoeff;

private:
  std::unique_ptr< constitutive::compositional::ComponentProperties > m_component_properties{};

private:
  template< typename ARRAY, typename LIST, typename DATAARRAY >
  static void createArray( ARRAY & array, LIST const & indices, integer const row, DATAARRAY const & data )
  {
    for( auto const i : indices )
    {
      array.emplace_back( data[Fluid::NC *row+i] );
    }
  }
public:
  template< typename ARRAY, typename LIST >
  static void createArray( ARRAY & array, LIST const & data )
  {
    for( auto const value : data )
    {
      array.emplace_back( value );
    }
  }
  template< typename ARRAY, typename LIST >
  static void populateArray( ARRAY & array, LIST const & data )
  {
    integer i = 0;
    for( auto const value : data )
    {
      array[i++] = value;
    }
  }
};

std::array< real64, 66 > Fluid::data = {
  // -- Pc
  2.2050e+07, 7.3750e+06, 3.4000e+06, 8.9630e+06, 1.2960e+06, 4.8721e+06,
  4.2481e+06, 3.6400e+06, 4.5990e+06, 2.5300e+06, 1.4600e+06,
  // -- Tc
  6.4700e+02, 3.0410e+02, 1.2620e+02, 3.7353e+02, 3.3150e+01, 3.0532e+02,
  3.6983e+02, 4.0785e+02, 1.9060e+02, 6.2200e+02, 7.8200e+02,
  // -- Vc
  6.4920e-05, 9.1025e-05, 8.1615e-05, 9.2053e-05, 5.5585e-05, 1.3810e-04,
  1.9170e-04, 2.4649e-04, 9.1302e-05, 5.3923e-04, 1.1664e-03,
  // -- Ac
  3.4400e-01, 2.3900e-01, 4.0000e-02, 9.4200e-02, -2.1900e-01, 9.9500e-02,
  1.5230e-01, 1.8440e-01, 1.1400e-02, 4.4300e-01, 8.1600e-01,
  // -- Mw
  1.8015e-02, 4.4010e-02, 2.8013e-02, 3.4100e-02, 1.6043e-02, 3.0070e-02,
  4.4097e-02, 5.8124e-02, 7.2151e-02, 1.1423e-01, 1.4228e-01,
  // -- Vs
  0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
  0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
};

}// testing

}// geos

#endif //GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUID_HPP_
