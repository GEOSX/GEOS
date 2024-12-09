/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file KValueFlashParameters.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHPARAMETERS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHPARAMETERS_HPP_

#include "ModelParameters.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

class TableFunction;

namespace constitutive
{

namespace compositional
{

template< integer NUM_PHASE >
class KValueFlashParameters : public ModelParameters
{
  static constexpr integer numPhases = NUM_PHASE;
public:
  KValueFlashParameters( std::unique_ptr< ModelParameters > parameters );
  ~KValueFlashParameters() override = default;

  static std::unique_ptr< ModelParameters > create( std::unique_ptr< ModelParameters > parameters );

  void createTables( string const & name,
                     string & pressureTableName,
                     string & temperatureTableName ) const;

  array1d< real64 > m_pressureCoordinates;
  array1d< real64 > m_temperatureCoordinates;
  string_array m_kValueTables;

  array1d< array1d< real64 > > m_pressureValues;
  array1d< array1d< real64 > > m_temperatureValues;
  array4d< real64 > m_kValueHyperCube;

  struct viewKeyStruct
  {
    static constexpr char const * pressureCoordinatesString() { return "pressureCoordinates"; }
    static constexpr char const * temperatureCoordinatesString() { return "temperatureCoordinates"; }
    static constexpr char const * kValueTablesString() { return "kValueTables"; }
  };

protected:
  void registerParametersImpl( MultiFluidBase * fluid ) override;
  void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties ) override;

private:
  static bool isIncreasing( arraySlice1d< real64 const > const & array );

  void generateHyperCube( integer const numComps );
  bool validateKValues( MultiFluidBase const * fluid ) const;
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHPARAMETERS_HPP_
