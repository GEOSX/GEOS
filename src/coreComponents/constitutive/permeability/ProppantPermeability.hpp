/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ProppantPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_PROPPANTPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_PROPPANTPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geos
{
namespace constitutive
{

class ProppantPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  ProppantPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                              arrayView3d< real64 > const & dPerm_dPressure,
                              arrayView4d< real64 > const & dPerm_dDispJump,
                              arrayView3d< real64 > const & permeabilityMultiplier,
                              real64 const & proppantPackPermeability )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dDispJump( dPerm_dDispJump ),
    m_permeabilityMultiplier( permeabilityMultiplier ),
    m_proppantPackPermeability( proppantPackPermeability )
  {}

  GEOS_HOST_DEVICE
  void compute( real64 const & oldHydraulicAperture,
                real64 const & newHydraulicAperture,
                real64 const & proppantPackVolumeFraction,
                arraySlice1d< real64 > const & permeability,
                arraySlice2d< real64 > const & dPerm_dDispJump,
                arraySlice1d< real64 > const & permeabilityMultiplier ) const
  {

    // permeability
    real64 const perm = ( oldHydraulicAperture*oldHydraulicAperture*oldHydraulicAperture +
                          oldHydraulicAperture*oldHydraulicAperture*newHydraulicAperture +
                          oldHydraulicAperture*newHydraulicAperture*newHydraulicAperture +
                          newHydraulicAperture*newHydraulicAperture*newHydraulicAperture ) / 48.0;

    real64 const dPerm  = ( oldHydraulicAperture*oldHydraulicAperture +
                            2.0*oldHydraulicAperture*newHydraulicAperture +
                            3.0*newHydraulicAperture*newHydraulicAperture ) / 48.0;

    real64 const squaredHydraulicAperture = newHydraulicAperture * newHydraulicAperture;

    // horizontal multiplier
    permeabilityMultiplier[0] = ( 1.0 - proppantPackVolumeFraction ) + 12.0 * proppantPackVolumeFraction * m_proppantPackPermeability / squaredHydraulicAperture;

    // vertical multiplier
    permeabilityMultiplier[1] = 1.0 / (1.0 - proppantPackVolumeFraction + proppantPackVolumeFraction * squaredHydraulicAperture / ( 12.0 * m_proppantPackPermeability ) );

    for( int dim=0; dim < 3; dim++ )
    {
      permeability[dim]        = perm;
      dPerm_dDispJump[dim][0]  = dPerm;
      dPerm_dDispJump[dim][1]  = 0.0;
      dPerm_dDispJump[dim][2]  = 0.0;

    }
  }

  GEOS_HOST_DEVICE
  virtual void updateFromApertureAndProppantVolumeFraction ( localIndex const k,
                                                             localIndex const q,
                                                             real64 const & oldHydraulicAperture,
                                                             real64 const & newHydraulicAperture,
                                                             real64 const & dHydraulicAperture_dNormalJump,
                                                             real64 const & proppantPackVolumeFraction ) const override final
  {
    GEOS_UNUSED_VAR( q, dHydraulicAperture_dNormalJump );

    compute( oldHydraulicAperture,
             newHydraulicAperture,
             proppantPackVolumeFraction,
             m_permeability[k][0],
             m_dPerm_dDispJump[k][0],
             m_permeabilityMultiplier[k][0] );
  }

private:

  arrayView4d< real64 > m_dPerm_dDispJump;

  arrayView3d< real64 > m_permeabilityMultiplier;

  real64 m_proppantPackPermeability;

};


class ProppantPermeability : public PermeabilityBase
{
public:

  ProppantPermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ProppantPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ProppantPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dDispJump,
                          m_permeabilityMultiplier,
                          m_proppantPackPermeability );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * maxProppantConcentrationString() { return "maxProppantConcentration"; }
    static constexpr char const * proppantDiameterString() { return "proppantDiameter"; }

    static constexpr char const * proppantPackPermeabilityString() { return "proppantPackPermeability"; }

  };

protected:

  virtual void postInputInitialization() override;

private:

  array4d< real64 > m_dPerm_dDispJump;

  array3d< real64 > m_permeabilityMultiplier;

  real64 m_proppantDiameter;

  real64 m_maxProppantConcentration;

  real64 m_proppantPackPermeability;

};

} /* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_PROPPANTPERMEABILITY_HPP_
