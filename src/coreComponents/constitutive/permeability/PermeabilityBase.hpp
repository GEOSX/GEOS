/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PermeabilityBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_PERMEABILITYBASE_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_PERMEABILITYBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geos
{
namespace constitutive
{

class PermeabilityBaseUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_permeability.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_permeability.size( 1 ); }

  GEOS_HOST_DEVICE
  virtual void updateFromPressureAndPorosity( localIndex const k,
                                              localIndex const q,
                                              real64 const & pressure,
                                              real64 const & porosity ) const
  {
    GEOS_UNUSED_VAR( k, q, pressure, porosity );
  }

  GEOS_HOST_DEVICE
  virtual void updateFromAperture( localIndex const k,
                                   localIndex const q,
                                   real64 const & oldHydraulicAperture,
                                   real64 const & newHydraulicAperture,
                                   real64 const & dHydraulicAperture_dNormalJump ) const
  {
    GEOS_UNUSED_VAR( k, q, oldHydraulicAperture, newHydraulicAperture, dHydraulicAperture_dNormalJump );
  }

  GEOS_HOST_DEVICE
  virtual void updateFromApertureAndShearDisplacement( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & oldHydraulicAperture,
                                                       real64 const & newHydraulicAperture,
                                                       real64 const & dHydraulicAperture_dNormalJump,
                                                       real64 const & pressure,
                                                       real64 const ( &dispJump )[3],
                                                       real64 const ( &traction )[3] ) const
  {
    GEOS_UNUSED_VAR( k, q, oldHydraulicAperture, newHydraulicAperture, dHydraulicAperture_dNormalJump, dispJump, traction, pressure );
  }

  GEOS_HOST_DEVICE
  virtual void updateFromApertureAndProppantVolumeFraction ( localIndex const k,
                                                             localIndex const q,
                                                             real64 const & oldHydraulicAperture,
                                                             real64 const & newHydraulicAperture,
                                                             real64 const & dHydraulicAperture_dNormalJump,
                                                             real64 const & proppantPackVolumeFraction ) const
  {
    GEOS_UNUSED_VAR( k, q, oldHydraulicAperture, newHydraulicAperture, dHydraulicAperture_dNormalJump, proppantPackVolumeFraction );
  }

protected:

  PermeabilityBaseUpdate( arrayView3d< real64 > const & permeability,
                          arrayView3d< real64 > const & dPerm_dPressure )
    : m_permeability( permeability ),
    m_dPerm_dPressure( dPerm_dPressure )
  {}

  arrayView3d< real64 > m_permeability;

  arrayView3d< real64 > m_dPerm_dPressure;
};

class PermeabilityBase : public ConstitutiveBase
{
public:

  PermeabilityBase( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @brief Const/non-mutable accessor for permeability.
   * @return Accessor
   */
  arrayView3d< real64 const > permeability() const { return m_permeability; }

  /**
   * @brief Const/non-mutable accessor for dPerm_dPressure.
   * @return Accessor
   */
  arrayView3d< real64 const > dPerm_dPressure() const { return m_dPerm_dPressure; }

  /**
   * @brief Utility function to scale the horizontal permeability (for instance, by net-to-gross)
   * @param[in] scalingFactors the vector of scaling factors (one value per cell) for the horizontal permeability
   */
  void scaleHorizontalPermeability( arrayView1d< real64 const > scalingFactors ) const;

  /**
   * @brief Initialize the permeability state
   */
  virtual void initializeState() const
  {}

protected:

  /// Vector of absolute permeability
  array3d< real64 > m_permeability;

  /// Vector of derivative of permeability wrt pressure
  array3d< real64 > m_dPerm_dPressure;
};

}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_PERMEABILITYBASE_HPP_
