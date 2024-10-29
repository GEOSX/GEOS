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
 * @file FrictionCoefficientSwapMPMEvent.hpp
 */

#ifndef GEOSX_FRICTIONCOEFFICIENTSWAP_MPMEVENT_HPP_
#define GEOSX_FRICTIONCOEFFICIENTSWAP_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

namespace geos
{

/**
 * @class FrictionCoefficientSwapMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class FrictionCoefficientSwapMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  FrictionCoefficientSwapMPMEvent( const string & name,
                Group * const parent );

  /// Destructor
  virtual ~FrictionCoefficientSwapMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "FrictionCoefficientSwap"; }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * frictionCoefficientString() { return "frictionCoefficient"; }
    static constexpr char const * frictionCoefficientTableString() { return "frictionCoefficientTable"; }

  } FrictionCoefficientSwapMPMEventViewKeys;
  /// @endcond

    real64 getFrictionCoefficient() const { return m_frictionCoefficient; }
    array2d< real64 > getFrictionCoefficientTable() const { return m_frictionCoefficientTable; }

protected:
  virtual void postInputInitialization() override final;

  // Event variables
  real64 m_frictionCoefficient;
  array2d< real64 > m_frictionCoefficientTable; 
};

} /* namespace geos */

#endif /* GEOSX_FRICTIONCOEFFICIENTSWAP_MPMEVENT_HPP_ */
