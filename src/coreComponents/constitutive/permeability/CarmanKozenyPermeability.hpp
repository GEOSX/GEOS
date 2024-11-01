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
 * @file CarmanKozenyPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_CARMANKOZENYPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_CARMANKOZENYPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geos
{
namespace constitutive
{

class CarmanKozenyPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  CarmanKozenyPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                  arrayView3d< real64 > const & dPerm_dPressure,
                                  arrayView3d< real64 > const & dPerm_dPorosity,
                                  real64 const particleDiameter,
                                  real64 const sphericity,
                                  R1Tensor const anisotropy )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dPorosity( dPerm_dPorosity ),
    m_particleDiameter( particleDiameter ),
    m_sphericity( sphericity ),
    m_anisotropy( anisotropy )
  {}

  GEOS_HOST_DEVICE
  void compute( real64 const & porosity,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dPorosity ) const;

  GEOS_HOST_DEVICE
  virtual void updateFromPressureAndPorosity( localIndex const k,
                                              localIndex const q,
                                              real64 const & pressure,
                                              real64 const & porosity ) const override
  {
    GEOS_UNUSED_VAR( pressure );

    compute( porosity,
             m_permeability[k][q],
             m_dPerm_dPorosity[k][q] );
  }

private:

  /// dPermeability_dPorosity
  arrayView3d< real64 > m_dPerm_dPorosity;

  /// Particle diameter
  real64 m_particleDiameter;

  /// Sphericity of the particles
  real64 m_sphericity;

  /// Anisotropy factors for three dimensions
  R1Tensor m_anisotropy;
};


class CarmanKozenyPermeability : public PermeabilityBase
{
public:

  CarmanKozenyPermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "CarmanKozenyPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CarmanKozenyPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dPorosity,
                          m_particleDiameter,
                          m_sphericity,
                          m_anisotropy );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * dPerm_dPorosityString() { return "dPerm_dPorosity"; }
    static constexpr char const * particleDiameterString() { return "particleDiameter"; }
    static constexpr char const * sphericityString() { return "sphericity"; }
    static constexpr char const * anisotropyString() { return "anisotropy"; }
  } viewKeys;

private:

  /// dPermeability_dPorosity
  array3d< real64 > m_dPerm_dPorosity;

  /// Particle diameter
  real64 m_particleDiameter;

  /// Sphericity of the particles
  real64 m_sphericity;

  /// Anisotropy factors for three dimensions
  R1Tensor m_anisotropy;
};


GEOS_HOST_DEVICE
inline
void CarmanKozenyPermeabilityUpdate::compute( real64 const & porosity,
                                              arraySlice1d< real64 > const & permeability,
                                              arraySlice1d< real64 > const & dPerm_dPorosity ) const
{
  real64 const constant = pow( m_sphericity*m_particleDiameter, 2 ) / 180;

  real64 const permValue = constant * pow( porosity, 3 )/ pow( (1 - porosity), 2 );

  real64 const dPerm_dPorValue = -constant * ( (porosity - 3) *  pow( porosity, 2 ) / pow( (1-porosity), 3 )  );

  for( localIndex i = 0; i < permeability.size(); ++i )
  {
    permeability[i] = permValue * m_anisotropy[i];
    dPerm_dPorosity[i] = dPerm_dPorValue * m_anisotropy[i];
  }
}



}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_CARMANKOZENYPERMEABILITY_HPP_
