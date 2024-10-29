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
 * @file CrystalHealMPMEvent.hpp
 */

#ifndef GEOSX_CRYSTALHEAL_MPMEVENT_HPP_
#define GEOSX_CRYSTALHEAL_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class CrystalHealMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class CrystalHealMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  CrystalHealMPMEvent( const string & name,
                       Group * const parent );

  /// Destructor
  virtual ~CrystalHealMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "CrystalHeal"; }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * targetRegionString() { return "targetRegion"; }
    static constexpr char const * healTypeString() { return "healType"; }
    static constexpr char const * markedParticlesToHealString() { return "markedParticlesToHeal"; }

    dataRepository::ViewKey targetRegion = { targetRegionString() };
  } CrystalHealMPMEventViewKeys;
  /// @endcond

  string getTargetRegion() const { return m_targetRegion; }
  int getHealType() const { return m_healType; }
  int getMarkedParticlesToHeal() const { return m_markedParticlesToHeal; }

  void setMarkedParticlesToHeal(int markedParticlesToHeal ) { m_markedParticlesToHeal = markedParticlesToHeal; }

protected:
  virtual void postInputInitialization() override final;

  // Event variables
  string m_targetRegion;
  int m_healType;
  int m_markedParticlesToHeal;
};

} /* namespace geos */

#endif /* GEOSX_CRYSTALHEAL_MPMEVENT_HPP_ */
