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
 * @file PressurePermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_PRESSUREPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_PRESSUREPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geos
{
namespace constitutive
{

enum class PressureModelType : integer
{
  Exponential,
  Hyperbolic
};

ENUM_STRINGS( PressureModelType,
              "Exponential",
              "Hyperbolic" );

class PressurePermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  PressurePermeabilityUpdate( PressureModelType const & presModelType,
                              R1Tensor const pressureDependenceConstants,
                              real64 const & referencePressure,
                              real64 const & maxPermeability,
                              arrayView3d< real64 > const & referencePermeability,
                              arrayView3d< real64 > const & permeability,
                              arrayView3d< real64 > const & dPerm_dPressure )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_presModelType( presModelType ),
    m_pressureDependenceConstants( pressureDependenceConstants ),
    m_referencePressure( referencePressure ),
    m_maxPermeability( maxPermeability ),
    m_referencePermeability( referencePermeability )
  {}

  GEOS_HOST_DEVICE
  void compute( real64 const & deltaPressure,
                R1Tensor const pressureDependenceConstants,
                real64 const (&referencePermeability)[3],
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dPressure ) const;

  GEOS_HOST_DEVICE
  void compute( real64 const & deltaPressure,
                R1Tensor const pressureDependenceConstants,
                real64 const (&referencePermeability)[3],
                real64 const maxPermeability,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dPressure ) const;

  GEOS_HOST_DEVICE
  virtual void updateFromPressureAndPorosity( localIndex const k,
                                              localIndex const q,
                                              real64 const & pressure,
                                              real64 const & porosity ) const override
  {
    GEOS_UNUSED_VAR( q, porosity );

    real64 const deltaPressure = pressure - m_referencePressure;

    real64 referencePermeability[3];

    referencePermeability[0] = m_referencePermeability[k][0][0];
    referencePermeability[1] = m_referencePermeability[k][0][1];
    referencePermeability[2] = m_referencePermeability[k][0][2];

    switch( m_presModelType )
    {
      case PressureModelType::Exponential:
      {
        compute( deltaPressure,
                 m_pressureDependenceConstants,
                 referencePermeability,
                 m_permeability[k][0],
                 m_dPerm_dPressure[k][0] );

        break;
      }
      case PressureModelType::Hyperbolic:
      {
        compute( deltaPressure,
                 m_pressureDependenceConstants,
                 referencePermeability,
                 m_maxPermeability,
                 m_permeability[k][0],
                 m_dPerm_dPressure[k][0] );

        break;
      }
      default:
      {
        GEOS_ERROR( "PressureModelType is invalid! It should be either Exponential or Hyperbolic" );
      }
    }
  }

private:

  /// Pressure dependence model type
  PressureModelType m_presModelType;

  /// Pressure dependent coefficients for each permeability component
  R1Tensor m_pressureDependenceConstants;

  /// Reference pressure in the model
  real64 const m_referencePressure;

  /// Maximum permeability
  real64 const m_maxPermeability;

  arrayView3d< real64 > m_referencePermeability;

};


class PressurePermeability : public PermeabilityBase
{
public:

  PressurePermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PressurePermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = PressurePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_presModelType,
                          m_pressureDependenceConstants,
                          m_referencePressure,
                          m_maxPermeability,
                          m_referencePermeability,
                          m_permeability,
                          m_dPerm_dPressure );
  }

  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * referencePermeabilityComponentsString() { return "referencePermeabilityComponents"; }
    static constexpr char const * pressureDependenceConstantsString() { return "pressureDependenceConstants"; }
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * referencePermeabilityString() { return "referencePermeability"; }
    static constexpr char const * maxPermeabilityString() { return "maxPermeability"; }
    static constexpr char const * pressureModelTypeString() { return "pressureModelType"; }
  } viewKeys;

  virtual void initializeState() const override final;

protected:

  virtual void postInputInitialization() override;

private:

  /// Permeability components at the reference pressure
  R1Tensor m_referencePermeabilityComponents;

  /// Pressure dependent coefficients for each permeability component
  R1Tensor m_pressureDependenceConstants;

  /// Reference pressure in the model
  real64 m_referencePressure;

  /// Maximum permeability
  real64 m_maxPermeability;

  array3d< real64 > m_referencePermeability;

  /// Pressure dependence model type
  PressureModelType m_presModelType;

};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void PressurePermeabilityUpdate::compute( real64 const & deltaPressure,
                                          R1Tensor const pressureDependenceConstants,
                                          real64 const (&referencePermeability)[3],
                                          arraySlice1d< real64 > const & permeability,
                                          arraySlice1d< real64 > const & dPerm_dPressure ) const
{
  for( localIndex i=0; i < permeability.size(); i++ )
  {
    real64 const perm = referencePermeability[i] * exp( pressureDependenceConstants[i] * deltaPressure );

    permeability[i] = perm;
    dPerm_dPressure[i] = perm * pressureDependenceConstants[i];
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void PressurePermeabilityUpdate::compute( real64 const & deltaPressure,
                                          R1Tensor const pressureDependenceConstants,
                                          real64 const (&referencePermeability)[3],
                                          real64 const maxPermeability,
                                          arraySlice1d< real64 > const & permeability,
                                          arraySlice1d< real64 > const & dPerm_dPressure ) const
{
  for( localIndex i=0; i < permeability.size(); i++ )
  {
    real64 const pressureOffSet = log( maxPermeability/referencePermeability[i] - 1 )/pressureDependenceConstants[i];

    real64 const perm = maxPermeability/( 1 + exp( -pressureDependenceConstants[i]*( deltaPressure - pressureOffSet ) ) );
    permeability[i] = perm;
    dPerm_dPressure[i] = perm*perm/maxPermeability*pressureDependenceConstants[i]*exp( -pressureDependenceConstants[i]*deltaPressure );
  }
}

}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_PRESSUREPERMEABILITY_HPP_
