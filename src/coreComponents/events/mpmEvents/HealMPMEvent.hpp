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
 * @file HealMPMEvent.hpp
 */

#ifndef GEOSX_HEAL_MPMEVENT_HPP_
#define GEOSX_HEAL_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class HealMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class HealMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  HealMPMEvent( const string & name,
                Group * const parent );

  /// Destructor
  virtual ~HealMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "Heal"; }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * targetRegionString() { return "targetRegion"; }

    dataRepository::ViewKey targetRegion = { targetRegionString() };

  } HealMPMEventViewKeys;
  /// @endcond

protected:
  virtual void postInputInitialization() override final;

  // Event variables
  string m_targetRegion;
};

} /* namespace geos */

#endif /* GEOSX_HEAL_MPMEVENT_HPP_ */
