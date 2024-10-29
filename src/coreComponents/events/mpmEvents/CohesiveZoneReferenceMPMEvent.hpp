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
 * @file CohesiveZoneReferenceMPMEvent.hpp
 */

#ifndef GEOSX_COHESIVEZONEREFERENCE_MPMEVENT_HPP_
#define GEOSX_COHESIVEZONEREFERENCE_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class CohesiveZoneReferenceMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class CohesiveZoneReferenceMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  CohesiveZoneReferenceMPMEvent( const string & name,
                  Group * const parent );

  /// Destructor
  virtual ~CohesiveZoneReferenceMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "CohesiveZoneReference"; }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {

  } CohesiveZoneReferenceMPMEventViewKeys;
  /// @endcond

protected:
  virtual void postInputInitialization() override final;

};

} /* namespace geos */

#endif /* GEOSX_COHESIVEZONEREFERENCE_MPMEVENT_HPP_ */
