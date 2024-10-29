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
 * @file BodyForceUpdateMPMEvent.hpp
 */

#ifndef GEOSX_BODYFORCEUPDATE_MPMEVENT_HPP_
#define GEOSX_BODYFORCEUPDATE_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class BodyForceUpdateMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class BodyForceUpdateMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  BodyForceUpdateMPMEvent( const string & name,
                  Group * const parent );

  /// Destructor
  virtual ~BodyForceUpdateMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "BodyForceUpdate"; }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * bodyForceString() { return "bodyForce"; }
  } BodyForceUpdateMPMEventViewKeys;
  /// @endcond

  array1d< real64 > getBodyForce() const { return m_bodyForce; } 

protected:
  virtual void postInputInitialization() override final;

  // Event variables
  array1d< real64 > m_bodyForce;
};

} /* namespace geos */

#endif /* GEOSX_BODYFORCEUPDATE_MPMEVENT_HPP_ */
