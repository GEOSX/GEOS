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
 * @file InsertPeriodicContactSurfacesMPMEvent.hpp
 */

#ifndef GEOSX_INSERTPERIODICCONTACTSURFACES_MPMEVENT_HPP_
#define GEOSX_INSERTPERIODICCONTACTSURFACES_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class InsertPeriodicContactSurfacesMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class InsertPeriodicContactSurfacesMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  InsertPeriodicContactSurfacesMPMEvent(  const string & name,
                         Group * const parent );

  /// Destructor
  virtual ~InsertPeriodicContactSurfacesMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "InsertPeriodicContactSurfaces"; }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {

  } InsertPeriodicContactSurfacesMPMEventViewKeys;
  /// @endcond

};

} /* namespace geos */

#endif /* GEOSX_INSERTPERIODICCONTACTSURFACES_MPMEVENT_HPP_ */
