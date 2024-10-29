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
 * @file MaterialSwapMPMEvent.hpp
 */

#ifndef GEOSX_MATERIALSWAP_MPMEVENT_HPP_
#define GEOSX_MATERIALSWAP_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class MaterialSwapMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class MaterialSwapMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  MaterialSwapMPMEvent(  const string & name,
                         Group * const parent );

  /// Destructor
  virtual ~MaterialSwapMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "MaterialSwap"; }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * sourceRegionString() { return "sourceRegion"; }
    static constexpr char const * destinationRegionString() { return "destinationRegion"; }

    dataRepository::ViewKey sourceRegion = { sourceRegionString() };
    dataRepository::ViewKey destinationRegion = { destinationRegionString() };

  } materialSwapMPMEventViewKeys;
  /// @endcond

  string getSourceRegion() const { return m_sourceRegion; }
  string getDestinationRegion() const { return m_destinationRegion; }

protected:
  virtual void postInputInitialization() override final;

  // Event variables
  string m_sourceRegion;
  string m_destinationRegion;
};

} /* namespace geos */

#endif /* GEOSX_MATERIALSWAP_MPMEVENT_HPP_ */
