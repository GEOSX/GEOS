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
 * @file EquilibriumInitialCondition.hpp
 */


#ifndef GEOS_FIELDSPECIFICATION_EQUILIBRIUMINITIALCONDITION_HPP
#define GEOS_FIELDSPECIFICATION_EQUILIBRIUMINITIALCONDITION_HPP

#include "FieldSpecificationBase.hpp"

namespace geos
{

/**
 * @class EquilibriumInitialCondition
 * Holds data to compute an hydrostatic equilibrium condition for flow problems
 */
class EquilibriumInitialCondition : public FieldSpecificationBase
{
public:

  /// @copydoc FieldSpecificationBase(string const &, dataRepository::Group *)
  EquilibriumInitialCondition( string const & name, Group * parent );

  /// deleted default constructor
  EquilibriumInitialCondition() = delete;

  /// default destructor
  virtual ~EquilibriumInitialCondition() = default;

  /// deleted copy constructor
  EquilibriumInitialCondition( EquilibriumInitialCondition const & ) = delete;

  /// defaulted move constructor
  EquilibriumInitialCondition( EquilibriumInitialCondition && ) = default;

  /// deleted copy assignment operator
  EquilibriumInitialCondition & operator=( EquilibriumInitialCondition const & ) = delete;

  /// deleted move assignment operator
  EquilibriumInitialCondition & operator=( EquilibriumInitialCondition && ) = delete;

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "HydrostaticEquilibrium"; }

  /**
   * @brief Getter for the max number of equilibration iterations
   * @return the max number of equilibrium iterations
   */
  integer getMaxNumEquilibrationIterations() const { return m_maxNumEquilibrationIterations; }

  /**
   * @brief Getter for the equilibration tolerance
   * @return the equilibration tolerance
   */
  real64 getEquilibrationTolerance() const { return m_equilibrationTolerance; }

  /**
   * @brief Getter for the elevation increment in the hydrostatic pressure table
   * @return the elevation increment in the hydrostatic pressure table
   */
  real64 getElevationIncrement() const { return m_elevationIncrement; }

  /**
   * @brief Getter for the datum pressure
   * @return the datum pressure
   */
  real64 getDatumPressure() const { return m_datumPressure; }

  /**
   * @brief Getter for the datum elevation
   * @return the datum elevation
   */
  real64 getDatumElevation() const { return m_datumElevation; }

  /**
   * @brief Getter for the component names
   * @return an array storing the component names
   */
  arrayView1d< string const > getComponentNames() const { return m_componentNames.toViewConst(); }

  /**
   * @brief Getter for the name of the phase initially saturating the reservoir
   * @return the name of the phase initially saturating the reservoir
   */
  string getInitPhaseName() const { return m_initPhaseName; }

  /**
   * @brief Getter for the component fraction table names
   * @return the component fraction table names
   */
  arrayView1d< string const > getComponentFractionVsElevationTableNames() const { return m_componentFractionVsElevationTableNames.toViewConst(); }

  /**
   * @brief Getter for the temperature table name
   * @return the temperature table name
   */
  string getTemperatureVsElevationTableName() const { return m_temperatureVsElevationTableName; }

  /**
   * @brief View keys
   */
  struct viewKeyStruct : public FieldSpecificationBase::viewKeyStruct
  {

    // equilibration parameters

    /// @return String key for the maximum number of equilibration iterations
    constexpr static char const * maxNumEquilibrationIterationsString() { return "maxNumberOfEquilibrationIterations"; }

    /// @return String key for the elevation increment in the hydrostatic pressure table
    constexpr static char const * elevationIncrementString() { return "elevationIncrementInHydrostaticPressureTable"; }

    /// @return String key for the equilibrium tolerance
    constexpr static char const * equilibrationToleranceString() { return "equilibrationTolerance"; }

    // datum elevation and pressure

    /// @return String key for the datum elevation
    constexpr static char const * datumElevationString() { return "datumElevation"; }

    /// @return String key for the datum pressure
    constexpr static char const * datumPressureString() { return "datumPressure"; }

    // name of the phase initially saturating the reservoir

    /// @return String key for the initial phase name
    constexpr static char const * initPhaseNameString() { return "initialPhaseName"; }

    // component names to use the component fraction tables

    /// @return String key for the component names
    constexpr static char const * componentNamesString() { return "componentNames"; }

    // tables storing the properties vs elevation properties

    /// @return String key for the component fraction vs elevation table names
    constexpr static char const * componentFractionVsElevationTableNamesString() { return "componentFractionVsElevationTableNames"; }

    /// @return String key for the temperature vs elevation table name
    constexpr static char const * temperatureVsElevationTableNameString() { return "temperatureVsElevationTableName"; }

  };


protected:

  virtual void postInputInitialization() override final;

  virtual void initializePreSubGroups() override final;

private:

  /// Maximum number of equilibration iterations
  integer m_maxNumEquilibrationIterations;

  /// Elevation increment in the hydrostatic pressure table
  real64 m_elevationIncrement;

  /// Tolerance used in the equilibration calculations
  real64 m_equilibrationTolerance;

  /// Datum (reference) elevation
  real64 m_datumElevation;

  /// Value of pressure at the datum elevation
  real64 m_datumPressure;

  /// Name of the phase initially saturating the reservoir
  string m_initPhaseName;

  /// Array of component names
  array1d< string > m_componentNames;

  /// Array of table names for component fraction vs elevation
  array1d< string > m_componentFractionVsElevationTableNames;

  /// Table name for temperature vs elevation
  string m_temperatureVsElevationTableName;

};

} /* namespace geos */

#endif /* GEOS_FIELDSPECIFICATION_EQUILIBRIUMINITIALCONDITION_HPP */
