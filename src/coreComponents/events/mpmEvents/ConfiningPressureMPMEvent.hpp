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
 * @file ConfiningPressureMPMEvent.hpp
 */

#ifndef GEOSX_CONFININGPRESSURE_MPMEVENT_HPP_
#define GEOSX_CONFININGPRESSURE_MPMEVENT_HPP_

#include "MPMEventBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

/**
 * @class ConfiningPressureMPMEvent
 *
 * This class implements a virtual fluid boundary condition for a region in the X-Y plane defined by a borehole radius.
 */
class ConfiningPressureMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  ConfiningPressureMPMEvent( const string & name,
                  Group * const parent );

  /// Destructor
  virtual ~ConfiningPressureMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "ConfiningPressure"; }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * confiningPressureBoxMinString() { return "confiningPressureBoxMin"; }
    static constexpr char const * confiningPressureBoxMaxString() { return "confiningPressureBoxMax"; }
    static constexpr char const * startPressureString() { return "startPressure"; }
    static constexpr char const * endPressureString() { return "endPressure"; }
    static constexpr char const * interpTypeString() { return "interpType"; }
  } ConfiningPressureMPMEventViewKeys;
  /// @endcond

  array1d< real64 > getConfiningPressureBoxMin() const { return m_confiningPressureBoxMin; } 
  array1d< real64 > getConfiningPressureBoxMax() const { return m_confiningPressureBoxMax; } 
  real64 getStartPressure() const { return m_startPressure; } 
  real64 getEndPressure() const { return m_endPressure; } 
  int getInterpType() const { return m_interpType; } 

protected:
  virtual void postInputInitialization() override final;

  // Event variables
  array1d< real64 > m_confiningPressureBoxMin;
  array1d< real64 > m_confiningPressureBoxMax;
  real64 m_startPressure;
  real64 m_endPressure;
  int m_interpType;

};

} /* namespace geos */

#endif /* GEOSX_CONFININGPRESSURE_MPMEVENT_HPP_ */
