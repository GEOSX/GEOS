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

/*
 * @file WellControls.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP

#include "codingUtilities/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{
namespace dataRepository
{
namespace keys
{
static constexpr auto wellControls = "WellControls";
}
}


/**
 * @class WellControls
 * @brief This class describes the controls used to operate a well.
 */
class WellControls : public dataRepository::Group
{
public:

  /** Type of wells
   * Either producer or injector.
   */
  enum class Type : integer
  {
    PRODUCER,  /**< A production well */
    INJECTOR   /**< An injection well */
  };


  /** Types of well controls
   * Used to specifiy a well's operating conditions
   */
  enum class Control : integer
  {
    BHP,  /**< The well operates at a specified bottom hole pressure (BHP) */
    PHASEVOLRATE, /**< The well operates at a specified phase volumetric flow rate */
    TOTALVOLRATE, /**< The well operates at a specified total volumetric flow rate */
    MASSRATE, /**<The well operates at a specified mass rate */
    UNINITIALIZED, /**< This is the current well control before postInputInitialization (needed to restart from file properly) */
  };


  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for WellControls Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit WellControls( string const & name, dataRepository::Group * const parent );


  /**
   * @brief Default destructor.
   */
  ~WellControls() override;

  /**
   * @brief Deleted default constructor.
   */
  WellControls() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  WellControls( WellControls const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  WellControls( WellControls && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a perforation object
   */
  WellControls & operator=( WellControls const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a perforation object
   */
  WellControls & operator=( WellControls && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Set the control type to BHP and set a numerical value for the control.
   * @param[in] val value for the BHP control
   */
  void switchToBHPControl( real64 const & val );

  /**
   * @brief Set the control type to total rate and set a numerical value for the control.
   * @param[in] val value for the total volumetric rate
   */
  void switchToTotalRateControl( real64 const & val );

  /**
   * @brief Set the control type to mass rate and set a numerical value for the control.
   * @param[in] val value for the mass rate
   */
  void switchToMassRateControl( real64 const & val );

  /**
   * @brief Set the control type to phase rate and set a numerical value for the control.
   * @param[in] val value for the phase volumetric rate
   */
  void switchToPhaseRateControl( real64 const & val );

  /**
   * @brief Get the control type for the well.
   * @return the Control enum enforced at the well
   */
  Control getControl() const { return m_currentControl; }

  /**
   * @brief Getter for the reference elevation where the BHP control is enforced
   * @return the reference elevation
   */
  real64 getReferenceElevation() const { return m_refElevation; }

  /**
   * @brief Getter for the reference gravity coefficient
   * @return the reference gravity coefficient
   */
  real64 getReferenceGravityCoef() const { return m_refGravCoef; }

  /**
   * @brief Setter for the reference gravity
   */
  void setReferenceGravityCoef( real64 const & refGravCoef ) { m_refGravCoef = refGravCoef; }


  /**
   * @brief Get the target bottom hole pressure value.
   * @return a value for the target bottom hole pressure
   */
  real64 getTargetBHP( real64 const & currentTime ) const
  {
    return m_targetBHPTable->evaluate( &currentTime );
  }

  /**
   * @brief Get the target total rate
   * @return the target total rate
   */
  real64 getTargetTotalRate( real64 const & currentTime ) const
  {
    return m_rateSign * m_targetTotalRateTable->evaluate( &currentTime );
  }

  /**
   * @brief Get the target phase rate
   * @return the target phase rate
   */
  real64 getTargetPhaseRate( real64 const & currentTime ) const
  {
    return m_rateSign * m_targetPhaseRateTable->evaluate( &currentTime );
  }

  /**
   * @brief Get the target mass rate
   * @return the target mass rate
   */
  real64 getTargetMassRate( real64 const & currentTime ) const
  {
    return m_rateSign * m_targetMassRateTable->evaluate( &currentTime );
  }

  /**
   * @brief Get the target phase name
   * @return the target phase name
   */
  const string & getTargetPhaseName() const { return m_targetPhaseName; }

  /**
   * @brief Const accessor for the composition of the injection stream
   * @return a global component fraction vector
   */
  arrayView1d< real64 const > getInjectionStream() const { return m_injectionStream; }

  /**
   * @brief Const accessor for the temperature of the injection stream
   * @return the temperature of the injection stream
   */
  real64 getInjectionTemperature() const { return m_injectionTemperature; }

  /**
   * @brief Getter for the flag specifying whether we check rates at surface or reservoir conditions
   * @return 1 if we use surface conditions, and 0 otherwise
   */
  integer useSurfaceConditions() const { return m_useSurfaceConditions; }

  /**
   * @brief Getter for the surface pressure when m_useSurfaceConditions == 1
   * @return the surface pressure
   */
  const real64 & getSurfacePressure() const { return m_surfacePres; }

  /**
   * @brief Getter for the surface temperature when m_useSurfaceConditions == 1
   * @return the surface temperature
   */
  const real64 & getSurfaceTemperature() const { return m_surfaceTemp; }

  /**
   * @brief Is the well an injector?
   * @return a boolean
   */
  bool isInjector() const { return ( m_type == Type::INJECTOR ); }

  /**
   * @brief Is the well a producer?
   * @return a boolean
   */
  bool isProducer() const { return ( m_type == Type::PRODUCER ); }

  /**
   * @brief Is the well open (or shut) at @p currentTime?
   * @param[in] currentTime the current time
   * @return a boolean
   */
  bool isWellOpen( real64 const & currentTime ) const;

  /**
   * @brief Getter for the flag to enable crossflow
   * @return the flag deciding whether crossflow is allowed or not
   */
  bool isCrossflowEnabled() const { return m_isCrossflowEnabled; }

  /**
   * @brief Getter for the initial pressure coefficient
   * @return the initial pressure coefficient
   */
  real64 getInitialPressureCoefficient() const { return m_initialPressureCoefficient; }

  /**
   * @brief set next time step based on tables intervals
   * @param[in] currentTime the current time
   * @param[inout] nextDt the time step
   */
  void setNextDtFromTables( real64 const currentTime, real64 & nextDt );

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the well reference elevation (for BHP control)
    static constexpr char const * refElevString() { return "referenceElevation"; }
    /// String key for the well type
    static constexpr char const * typeString() { return "type"; }
    /// String key for the well input control
    static constexpr char const * inputControlString() { return "control"; }
    /// String key for the well current control
    static constexpr char const * currentControlString() { return "currentControl"; }
    /// String key for the well target BHP
    static constexpr char const * targetBHPString() { return "targetBHP"; }
    /// String key for the well target rate
    static constexpr char const * targetTotalRateString() { return "targetTotalRate"; }
    /// String key for the well target phase rate
    static constexpr char const * targetPhaseRateString() { return "targetPhaseRate"; }
    /// String key for the well target phase name
    static constexpr char const * targetPhaseNameString() { return "targetPhaseName"; }
    /// String key for the well target phase name
    static constexpr char const * targetMassRateString() { return "targetMassRate"; }
    /// String key for the well injection stream
    static constexpr char const * injectionStreamString() { return "injectionStream"; }
    /// String key for the well injection temperature
    static constexpr char const * injectionTemperatureString() { return "injectionTemperature"; }
    /// String key for checking the rates at surface conditions
    static constexpr char const * useSurfaceConditionsString() { return "useSurfaceConditions"; }
    /// String key for the surface pressure
    static constexpr char const * surfacePressureString() { return "surfacePressure"; }
    /// String key for the surface temperature
    static constexpr char const * surfaceTemperatureString() { return "surfaceTemperature"; }
    /// string key for total rate table name
    static constexpr char const * targetTotalRateTableNameString() { return "targetTotalRateTableName"; }
    /// string key for phase rate table name
    static constexpr char const * targetPhaseRateTableNameString() { return "targetPhaseRateTableName"; }
    /// string key for mass rate table name
    static constexpr char const * targetMassRateTableNameString() { return "targetMassRateTableName"; }
    /// string key for BHP table name
    static constexpr char const * targetBHPTableNameString() { return "targetBHPTableName"; }
    /// string key for status table name
    static constexpr char const * statusTableNameString() { return "statusTableName"; }
    /// string key for the crossflow flag
    static constexpr char const * enableCrossflowString() { return "enableCrossflow"; }
    /// string key for the initial pressure coefficient
    static constexpr char const * initialPressureCoefficientString() { return "initialPressureCoefficient"; }

  }
  /// ViewKey struct for the WellControls class
  viewKeysWellControls;

protected:

  virtual void postInputInitialization() override;

  void setNextDtFromTable( TableFunction const * table, real64 const currentTime, real64 & nextDt );

private:

  /// Well type (as Type enum)
  Type m_type;

  /// Reference elevation
  real64 m_refElevation;

  /// Gravity coefficient of the reference elevation
  real64 m_refGravCoef;

  /// Input well controls as a Control enum
  Control m_inputControl;

  /// Well controls as a Control enum
  Control m_currentControl;

  /// Target bottom hole pressure value
  real64 m_targetBHP;

  /// Target rate value
  real64 m_targetTotalRate;

  /// Target phase rate value
  real64 m_targetPhaseRate;

  /// Name of the targeted phase
  string m_targetPhaseName;

  /// Target MassRate
  real64 m_targetMassRate;

  /// Vector with global component fractions at the injector
  array1d< real64 > m_injectionStream;

  /// Temperature at the injector
  real64 m_injectionTemperature;

  /// Flag to decide whether rates are controlled at rates or surface conditions
  integer m_useSurfaceConditions;

  /// Surface pressure
  real64 m_surfacePres;

  /// Surface temperature
  real64 m_surfaceTemp;

  /// Total rate table name
  string m_targetTotalRateTableName;

  /// Phase rate table name
  string m_targetPhaseRateTableName;

  /// Mass rate table name
  string m_targetMassRateTableName;

  /// BHP table name
  string m_targetBHPTableName;

  /// Status table name
  string m_statusTableName;

  /// Flag to enable crossflow
  integer m_isCrossflowEnabled;

  /// Tuning coefficient for the initial well pressure
  real64 m_initialPressureCoefficient;

  /// Rate sign. +1 for injector, -1 for producer
  real64 m_rateSign;

  /// Total rate table
  TableFunction const * m_targetTotalRateTable;

  /// Phase rate table
  TableFunction const * m_targetPhaseRateTable;

  /// Mass rate table
  TableFunction const * m_targetMassRateTable;

  /// BHP table
  TableFunction const * m_targetBHPTable;

  /// Status table
  TableFunction const * m_statusTable;
};

ENUM_STRINGS( WellControls::Type,
              "producer",
              "injector" );

ENUM_STRINGS( WellControls::Control,
              "BHP",
              "phaseVolRate",
              "totalVolRate",
              "massRate",
              "uninitialized" );


} //namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP
