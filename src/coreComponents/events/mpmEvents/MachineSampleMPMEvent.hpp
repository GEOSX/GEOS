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
 * @file MachineSampleMPMEvent.hpp
 */

#ifndef GEOSX_MACHINESAMPLE_MPMEVENT_HPP_
#define GEOSX_MACHINESAMPLE_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class MachineSampleMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class MachineSampleMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  MachineSampleMPMEvent(  const string & name,
                         Group * const parent );

  /// Destructor
  virtual ~MachineSampleMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "MachineSample"; }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * sampleTypeString() { return "sampleType"; }
    static constexpr char const * filletRadiusString() { return "filletRadius"; }
    static constexpr char const * gaugeLengthString() { return "gaugeLength"; }
    static constexpr char const * gaugeRadiusString() { return "gaugeRadius"; }
    static constexpr char const * diskRadiusString() { return "diskRadius"; }

    // dataRepository::ViewKey source = { sourceString() };
    // dataRepository::ViewKey destination = { destinationString() };

  } MachineSampleMPMEventViewKeys;
  /// @endcond

string getSampleType() const { return m_sampleType; }
real64 getFilletRadius() const { return m_filletRadius; }
real64 getGaugeLength() const { return m_gaugeLength; }
real64 getGaugeRadius() const { return m_gaugeRadius; }
real64 getDiskRadius() const { return m_diskRadius; }

protected:
  virtual void postInputInitialization() override final;

  // Event variables
  string m_sampleType;

  // Parameters for dogbone
  real64 m_filletRadius;
  real64 m_gaugeLength;
  real64 m_gaugeRadius;

  // Parameters for brazilDisk
  real64 m_diskRadius;
};

} /* namespace geos */

#endif /* GEOSX_MACHINESAMPLE_MPMEVENT_HPP_ */
