/*
 * -------------------------------------------------------------------------------------
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
 * @file Units.hpp
 * @brief Enumerates the Units that are in use in GEOS and regroups useful conversion and formatting functions.
 */
#ifndef GEOS_MATH_UNITS_HPP_
#define GEOS_MATH_UNITS_HPP_

#include "common/DataTypes.hpp"
#include "common/PhysicsConstants.hpp"

namespace geos
{

namespace units
{


/**
 * @return the input Kelvin degrees converted in Celsius
 * @param kelvin degrees input
 */
GEOS_HOST_DEVICE
inline constexpr double convertKToC( double kelvin )
{ return kelvin - constants::zeroDegreesCelsiusInKelvin; }
/**
 * @return the input Celsius degrees converted in Kelvin
 * @param celsius degrees input
 */
GEOS_HOST_DEVICE
inline constexpr double convertCToK( double celsius )
{ return celsius + constants::zeroDegreesCelsiusInKelvin; }


/**
 * @brief Enumerator of available unit types. Units are in SI by default.
 */
enum Unit : integer
{
  /// Default label when a unit is not known for a value
  Unknown,

  /// Label to use when a value has not physical dimension (ratio values, propotions...)
  Dimensionless,

  /// Pressure in Pascal
  Pressure,

  /// Temperature in Kelvin
  Temperature,

  /// Temperature in Celcius
  TemperatureInC,

  /// Distance in meter
  Distance,

  /// Time in seconds
  Time,

  /// Viscosity in Pa*s
  Viscosity,

  /// Density in kg/m³
  Density,

  /// Enthalpy in J/kg
  Enthalpy,

  /// Solubility in g/L
  Solubility,

  /// Mass in kg
  Mass,

  /// Mole in mol
  Mole,

  /// Mass rate in kg/s
  MassRate,

  /// Mole rate in mol/s
  MoleRate,
};


/**
 * @param unit The unit we want the information.
 * @return The name of the specified unit followed by its symbol in brackets.
 */
constexpr inline std::string_view getDescription( Unit unit )
{
  switch( unit )
  {
    default:              return "unknown [?]";
    case Dimensionless:   return "dimensionless [1]";
    case Pressure:        return "pressure [Pa]";
    case Temperature:     return "temperature [K]";
    case TemperatureInC:  return "temperature [C]";
    case Distance:        return "distance [m]";
    case Time:            return "time [s]";
    case Viscosity:       return "viscosity [Pa*s]";
    case Enthalpy:        return "enthalpy [J/kg]";
    case Density:         return "density [kg/m3]";
    case Solubility:      return "solubility [g/L]";
    case Mass:            return "mass [kg]";
    case Mole:            return "mole [mol]";
    case MassRate:        return "mass rate [kg/s]";
    case MoleRate:        return "mole rate [mol/s]";
  }
}

/**
 * @param unit The unit we want the information.
 * @return The symbol of the specified unit.
 */
constexpr inline std::string_view getSymbol( Unit unit )
{
  switch( unit )
  {
    default:              return "?";
    case Dimensionless:   return "1";
    case Pressure:        return "Pa";
    case Temperature:     return "K";
    case TemperatureInC:  return "C";
    case Distance:        return "m";
    case Time:            return "s";
    case Viscosity:       return "Pa*s";
    case Enthalpy:        return "J/kg";
    case Density:         return "kg/m3";
    case Solubility:      return "g/L";
    case Mass:            return "kg";
    case Mole:            return "mol";
    case MassRate:        return "kg/s";
    case MoleRate:        return "mol/s";
  }
}

/**
 * @brief Format the specified value coherently with the specified unit.
 * @param value The value to format.
 * @param unit The unit of the specified value.
 * @return A string that can be easily be integrated in a sentence and takes the form
 * "characteristic" of "value" ["unit symbol"]
 */
inline string formatValue( real64 value, Unit unit )
{
  switch( unit )
  {
    default:              return GEOS_FMT( "value of {} [?]", value );
    case Dimensionless:   return GEOS_FMT( "value of {} [1]", value );
    case Pressure:        return GEOS_FMT( "pressure of {} [Pa]", value );
    case Temperature:     return GEOS_FMT( "temperature of {} [K]", value );
    case TemperatureInC:  return GEOS_FMT( "temperature of {} [K]", convertCToK( value ) );
    case Distance:        return GEOS_FMT( "distance of {} [s]", value );
    case Time:            return GEOS_FMT( "time of {} [s]", value );
    case Viscosity:       return GEOS_FMT( "viscosity of {} [Pa*s]", value );
    case Enthalpy:        return GEOS_FMT( "enthalpy of {} [J/kg]", value );
    case Density:         return GEOS_FMT( "density of {} [kg/m3]", value );
    case Solubility:      return GEOS_FMT( "solubility of {} [g/L]", value );
    case Mass:            return GEOS_FMT( "mass of {} [kg]", value );
    case Mole:            return GEOS_FMT( "mole of {} [mol]", value );
    case MassRate:        return GEOS_FMT( "mass rate of {} [kg/s]", value );
    case MoleRate:        return GEOS_FMT( "mole rate of {} [mol/s]", value );
  }
}


/// Clock in use in GEOS to manipulate system times.
using SystemClock = std::chrono::system_clock;

/// One year = 365.2425 days (= 146097 / 400) following the Gregorian calendar and the C++ convention.
using YearDaysRatio = std::ratio< 146097, 400 >;
/// Day helper duration type, equivalent to C++20 std::chrono::days.
using Days = std::chrono::duration< int64_t, std::ratio_multiply< std::ratio< 24 >, std::chrono::hours::period > >;
/// Year helper duration type, equivalent to C++20 std::chrono::years.
using Years = std::chrono::duration< int64_t, std::ratio_multiply< YearDaysRatio, Days::period > >;

/// Days in one year (following the Gregorian calendar and the C++ convention) = 365.2425 days (= 146097 / 400).
static constexpr double YearDays = ( double )YearDaysRatio::num / YearDaysRatio::den;
/// Seconds in a minute
static constexpr double MinuteSeconds = 60.0;
/// Seconds in a hour
static constexpr double HourSeconds = 60.0 * MinuteSeconds;
/// Seconds in a day
static constexpr double DaySeconds = 24.0 * HourSeconds;
/// Seconds in a year
static constexpr double YearSeconds = YearDays * DaySeconds;


/**
 * @brief Stores information that is useful to duration strings. Based on the geos::units time constants
 */
struct TimeFormatInfo
{
  /// Total time (including the decimal part) this instance represents in seconds
  double const m_totalSeconds = 0.0;
  /// Number of integral years to show
  int const m_years = 0;
  /// Number of integral days to show
  int const m_days = 0;
  /// Number of integral hours to show
  int const m_hours = 0;
  /// Number of integral minutes to show
  int const m_minutes = 0;
  /// Number of integral seconds to show
  int const m_seconds = 0;

  /**
   * @brief Construct a TimeFormatInfo from raw data (which must be coherent)
   * @param totalSeconds The total time (including the decimal part) this instance represents in seconds
   * @param years Number of integral years to show
   * @param days Number of integral days to show
   * @param hours Number of integral hours to show
   * @param minutes Number of integral minutes to show
   * @param seconds Number of integral seconds to show
   */
  TimeFormatInfo( double totalSeconds, int years, int days, int hours, int minutes, int seconds );
  /**
   * @return A TimeFormatInfo constructed from the seconds to represent
   * @param seconds the total time to represents in seconds (including the decimal part)
   */
  static TimeFormatInfo fromSeconds( double const seconds );
  /**
   * @return A TimeFormatInfo constructed from a standard typed duration value
   * @param duration the duration to represents, in SystemClock::duration type
   * (more types can be added by adding std::chrono::duration template specialisations).
   */
  template< typename DURATION > static TimeFormatInfo fromDuration( DURATION duration );

  /**
   * @brief Insert the string representation information in the provided stream.
   */
  friend std::ostream & operator<<( std::ostream & os, TimeFormatInfo const & ctx );

  /**
   * @return a comprehensive string representation of this structure.
   */
  string toString() const;

  /**
   * @return a user friendly string representation of this structure (YDHMS format, without subsecond details).
   */
  string toUnfoldedString() const;

  /**
   * @return a precise string representation of this structure (decimal seconds).
   */
  string toSecondsString() const;
};


} // end namespace units

} // end namespace geos


/**
 * @brief Formatter to be able to directly use a DurationInfo as a GEOS_FMT() argument
 */
template<>
struct GEOS_FMT_NS::formatter< geos::units::TimeFormatInfo > : GEOS_FMT_NS::formatter< std::string >
{
  /**
   * @brief Format the specified TimeFormatInfo to a string.
   * @param durationData the TimeFormatInfo object to format
   * @param ctx formatting state consisting of the formatting arguments and the output iterator
   * @return iterator to the output buffer
   */
  auto format( geos::units::TimeFormatInfo const & durationData, format_context & ctx ) const
  {
    return GEOS_FMT_NS::formatter< std::string >::format( durationData.toString(), ctx );
  }
};

#endif //GEOS_MATH_PHYSICSCONSTANTS_HPP_
