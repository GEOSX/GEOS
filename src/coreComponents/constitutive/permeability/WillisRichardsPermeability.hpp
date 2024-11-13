/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file WillisRichardsPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_WILLISRICHARDSPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_WILLISRICHARDSPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geos
{
namespace constitutive
{

class WillisRichardsPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  WillisRichardsPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                    arrayView3d< real64 > const & dPerm_dPressure,
                                    arrayView4d< real64 > const & dPerm_dDispJump,
                                    arrayView4d< real64 > const & dPerm_dTraction,
                                    real64 const maxFracAperture,
                                    real64 const dilationCoefficient,
                                    real64 const refClosureStress )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dDispJump( dPerm_dDispJump ),
    m_dPerm_dTraction( dPerm_dTraction ),
    m_maxFracAperture( maxFracAperture ),
    m_dilationCoefficient( dilationCoefficient ),
    m_refClosureStress( refClosureStress )
  {}

  GEOS_HOST_DEVICE
  void compute( real64 const ( &dispJump )[3],
                real64 const ( &traction )[3],
                arraySlice1d< real64 > const & permeability,
                arraySlice2d< real64 > const & dPerm_dDispJump,
                arraySlice2d< real64 > const & dPerm_dTraction ) const;

  GEOS_HOST_DEVICE
  virtual void updateFromApertureAndShearDisplacement( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & oldHydraulicAperture,
                                                       real64 const & newHydraulicAperture,
                                                       real64 const & dHydraulicAperture_dNormalJump,
                                                       real64 const & pressure,
                                                       real64 const ( &dispJump )[3],
                                                       real64 const ( &traction )[3] ) const override final
  {
    GEOS_UNUSED_VAR( q, oldHydraulicAperture, newHydraulicAperture, dHydraulicAperture_dNormalJump, pressure );

    compute( dispJump,
             traction,
             m_permeability[k][0],
             m_dPerm_dDispJump[k][0],
             m_dPerm_dTraction[k][0] );
  }

private:

  /// Derivative of fracture permeability to shear displacement jump between fracture surfaces
  arrayView4d< real64 > m_dPerm_dDispJump;

  /// Derivative of fracture permeability to traction acting on fracture surfaces
  arrayView4d< real64 > m_dPerm_dTraction;

  /// Maximum fracture aperture at zero contact stress
  real64 m_maxFracAperture;

  /// Dilation coefficient
  real64 m_dilationCoefficient;

  /// Effective normal stress causes 90% reduction in aperture
  real64 m_refClosureStress;

};


class WillisRichardsPermeability : public PermeabilityBase
{
public:

  WillisRichardsPermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "WillisRichardsPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = WillisRichardsPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dDispJump,
                          m_dPerm_dTraction,
                          m_maxFracAperture,
                          m_dilationCoefficient,
                          m_refClosureStress );
  }

private:

  /// Derivative of fracture permeability to shear displacement jump between fracture surfaces
  array4d< real64 > m_dPerm_dDispJump;

  /// Derivative of fracture permeability to traction acting on fracture surfaces
  array4d< real64 > m_dPerm_dTraction;

  /// Maximum fracture aperture at zero contact stress
  real64 m_maxFracAperture;

  /// Dilation coefficient
  real64 m_dilationCoefficient;

  /// Effective normal stress causes 90% reduction in aperture
  real64 m_refClosureStress;

  struct viewKeyStruct
  {
    static constexpr char const * maxFracApertureString() { return "maxFracAperture"; }
    static constexpr char const * dilationCoefficientString() { return "dilationCoefficient"; }
    static constexpr char const * refClosureStressString() { return "refClosureStress"; }
  };

};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void WillisRichardsPermeabilityUpdate::compute( real64 const ( &dispJump )[3],
                                                real64 const ( &traction )[3],
                                                arraySlice1d< real64 > const & permeability,
                                                arraySlice2d< real64 > const & dPerm_dDispJump,
                                                arraySlice2d< real64 > const & dPerm_dTraction ) const
{
  real64 const shearMag = std::sqrt( dispJump[1]*dispJump[1] + dispJump[2]*dispJump[2] );

  real64 const effNormalStress = -traction[0];

  real64 const aperture = ( m_maxFracAperture + shearMag * m_dilationCoefficient ) / ( 1.0 + 9.0 * effNormalStress/m_refClosureStress );

  real64 const dPerm_daperture = aperture / 6.0;

  real64 const daperture_dshearMag = m_dilationCoefficient / ( 1.0 + 9.0 * effNormalStress/m_refClosureStress );

  real64 const daperture_deffNormalStress = -( m_maxFracAperture + shearMag * m_dilationCoefficient ) / ( 1.0 + 9.0 * effNormalStress/m_refClosureStress ) /
                                            ( 1.0 + 9.0 * effNormalStress/m_refClosureStress ) * 9.0 /m_refClosureStress;

  for( localIndex i=0; i < permeability.size(); i++ )
  {
    permeability[i] = aperture * aperture / 12.0;

    dPerm_dDispJump[i][0] = 0.0;
    real64 const tmpValue = shearMag > 0.0 ? daperture_dshearMag /shearMag : 0.0;
    dPerm_dDispJump[i][1] = dPerm_daperture * tmpValue * dispJump[1];
    dPerm_dDispJump[i][2] = dPerm_daperture * tmpValue * dispJump[2];

    dPerm_dTraction[i][0] = -dPerm_daperture * daperture_deffNormalStress;
    dPerm_dTraction[i][1] = 0.0;
    dPerm_dTraction[i][2] = 0.0;
  }
}



} /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
