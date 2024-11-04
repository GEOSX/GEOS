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
 * @file OBLFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_OBLFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_OBLFLUID_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "functions/MultivariableTableFunction.hpp"
#include "functions/PythonFunction.hpp"

namespace geos
{
namespace constitutive
{

enum OBLInterpolatorMode : integer
{
  Static = 0,         ///< static interpolation from a given tabled data
  Adaptive = 1        ///< adaptive interpolation from a given interface to function
};

enum OBLInterpolatorType : integer
{
  Multilinear = 0,    ///< multilinear interpolation
  Linear = 1          ///< linear interpolation
};

class OBLFluid : public ConstitutiveBase
{
public:
  typedef __uint128_t longIndex;

  OBLFluid( string const & name, Group * const parent );
  /**
   * @brief name of the constitutive in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string catalogName() { return "OBLFluid"; }
  /**
   * @copydoc ConstitutiveBase::getCatalogName()
   */
  virtual string getCatalogName() const override { return catalogName(); }
  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    // input
    static constexpr char const * interpolatorModeString() { return "interpolatorMode"; }
    static constexpr char const * interpolatorTypeString() { return "interpolatorType"; }
    static constexpr char const * oblOperatorsTableFileString() { return "oblOperatorsTableFile"; }
  };
  /**
   * @brief getter to the pointer to OBL operators table
   * @return pointer to OBL operators table.
   */
  MultivariableTableFunction const & getTable() const
  {
    GEOS_ERROR_IF( m_OBLOperatorsTable == nullptr, "m_OBLOperatorsTable is not initialized" );
    return *m_OBLOperatorsTable;
  }
  /**
   * @brief getter to the Python-based evaluator
   * @return pointer to the Python-based evaluator.
   */
  template< typename INDEX_T = longIndex >
  PythonFunction< INDEX_T > * getPythonFunction()
  {
    GEOS_ERROR_IF( m_pythonFunction == nullptr, "m_pythonFunction is not initialized" );
    return m_pythonFunction;
  }
  /**
   * @brief initialize input
   */
  void initialize()
  {
    if( !m_isInitialized )
    {
      postInputInitialization();
    }
  }
  /**
   * @brief Retrieves the current OBL interpolator mode.
   * @return OBLInterpolatorMode The current interpolation mode (Static or Adaptive).
   */
  OBLInterpolatorMode getInterpolatorMode() const { return m_interpolatorMode; };
  /**
   * @brief Retrieves the current OBL interpolator type.
   * @return OBLInterpolatorType The current interpolation type (Multilinear or Linear).
   */
  OBLInterpolatorType getInterpolatorType() const { return m_interpolatorType; };
private:
  /// OBL interpolator mode
  OBLInterpolatorMode m_interpolatorMode;

  /// corresponding input string
  string m_interpolatorModeString;

  /// OBL interpolator type
  OBLInterpolatorType m_interpolatorType;

  /// corresponding input string
  string m_interpolatorTypeString;

  /// OBL operators table file (if OBL physics becomes consitutive, multiple regions will be supported )
  Path m_OBLOperatorsTableFile;

  /// OBL operators table function tabulated vs all primary variables
  MultivariableTableFunction const * m_OBLOperatorsTable;

  /// OBL operators with access to Python-base exact evaluator
  PythonFunction< longIndex > * m_pythonFunction;

  /// Flag to check if contitutive is initialized or not
  bool m_isInitialized = false;

  /**
   * @copydoc dataRepository::Group::postInputInitialization()
   */
  virtual void postInputInitialization() override;
};

} //namespace constitutive

} //namespace geos


#endif //GEOS_CONSTITUTIVE_FLUID_OBLFLUID_HPP_
