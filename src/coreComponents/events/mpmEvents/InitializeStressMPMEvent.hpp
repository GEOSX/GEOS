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
 * @file InitializeStressMPMEvent.hpp
 */

#ifndef GEOSX_INITIALIZESTRESS_MPMEVENT_HPP_
#define GEOSX_INITIALIZESTRESS_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class InitializeStressMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class InitializeStressMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  InitializeStressMPMEvent( const string & name,
                  Group * const parent );

  /// Destructor
  virtual ~InitializeStressMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "InitializeStress"; }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * pressureString() { return "pressure"; }
    static constexpr char const * targetRegionString() { return "targetRegion"; }
   
    dataRepository::ViewKey targetRegion = { targetRegionString() };
  } InitializeStressMPMEventViewKeys;
  /// @endcond

  string getTargetRegion() const { return m_targetRegion; } 
  real64 getPressure() const { return m_pressure; } 

protected:
  virtual void postInputInitialization() override final;

  // Event variables
  real64 m_pressure;
  string m_targetRegion;

};

} /* namespace geos */

#endif /* GEOSX_INITIALIZESTRESS_MPMEVENT_HPP_ */
