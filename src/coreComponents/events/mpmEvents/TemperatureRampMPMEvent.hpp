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
 * @file TemperatureRampMPMEvent.hpp
 */

#ifndef GEOSX_TEMPERATURERAMP_MPMEVENT_HPP_
#define GEOSX_TEMPERATURERAMP_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class TemperatureRampMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class TemperatureRampMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  TemperatureRampMPMEvent( const string & name,
                  Group * const parent );

  /// Destructor
  virtual ~TemperatureRampMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "TemperatureRamp"; }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * startTemperatureString() { return "startTemperature"; }
    static constexpr char const * endTemperatureString() { return "endTemperature"; }
    static constexpr char const * interpTypeString() { return "interpType"; }
  } TemperatureRampMPMEventViewKeys;
  /// @endcond

  real64 getStartTemperature() const { return m_startTemperature; } 
  real64 getEndTemperature() const { return m_endTemperature; } 
  int getInterpType() const { return m_interpType; } 

protected:
  virtual void postInputInitialization() override final;

  // Event variables
  real64 m_startTemperature;
  real64 m_endTemperature;
  int m_interpType;

};

} /* namespace geos */

#endif /* GEOSX_TEMPERATURERAMP_MPMEVENT_HPP_ */
