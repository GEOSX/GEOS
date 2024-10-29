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
 * @file DeformationUpdateMPMEvent.hpp
 */

#ifndef GEOSX_DEFORMATIONUPDATE_MPMEVENT_HPP_
#define GEOSX_DEFORMATIONUPDATE_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class DeformationUpdateMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class DeformationUpdateMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  DeformationUpdateMPMEvent( const string & name,
                  Group * const parent );

  /// Destructor
  virtual ~DeformationUpdateMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "DeformationUpdate"; }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * prescribedFTableString() { return "prescribedFTable"; }
    static constexpr char const * prescribedBoundaryFTableString() { return "prescribedFTable"; }
    static constexpr char const * stressControlString() { return "stressControl"; }
  } DeformationUpdateMPMEventViewKeys;
  /// @endcond

  int getPrescribedBoundaryFTable() const { return m_prescribedBoundaryFTable; } 
  int getPrescribedFTable() const { return m_prescribedFTable; } 
  array1d< int > getStressControl() const { return m_stressControl; } 

protected:
  virtual void postInputInitialization() override final;

  // Event variables
  int m_prescribedFTable;
  int m_prescribedBoundaryFTable;
  array1d< int > m_stressControl;
};

} /* namespace geos */

#endif /* GEOSX_DEFORMATIONUPDATE_MPMEVENT_HPP_ */
