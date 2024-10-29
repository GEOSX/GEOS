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
 * @file ParallelPlatesPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geos
{
namespace constitutive
{

class ParallelPlatesPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  ParallelPlatesPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                    arrayView3d< real64 > const & dPerm_dPressure,
                                    arrayView4d< real64 > const & dPerm_dDispJump,
                                    bool const updateTransversalComponent )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dDispJump( dPerm_dDispJump ),
    m_numDimensionsToUpdate( 3 )
  {
    m_numDimensionsToUpdate = updateTransversalComponent ? 3 : 2;
  }

  GEOS_HOST_DEVICE
  void compute( real64 const & oldHydraulicAperture,
                real64 const & newHydraulicAperture,
                real64 const & dHydraulicAperture_dNormalJump,
                arraySlice1d< real64 > const & permeability,
                arraySlice2d< real64 > const & dPerm_dDispJump ) const
  {
    GEOS_UNUSED_VAR( oldHydraulicAperture );

    real64 const perm  = newHydraulicAperture*newHydraulicAperture / 12.0;
    real64 const dPerm_dHydraulicAperture = newHydraulicAperture / 6.0;

    for( int dim=0; dim < m_numDimensionsToUpdate; dim++ )
    {
      permeability[dim]        = perm;
      dPerm_dDispJump[dim][0]  = dPerm_dHydraulicAperture * dHydraulicAperture_dNormalJump;
      dPerm_dDispJump[dim][1]  = 0.0;
      dPerm_dDispJump[dim][2]  = 0.0;
    }
  }

  GEOS_HOST_DEVICE
  virtual void updateFromAperture( localIndex const k,
                                   localIndex const q,
                                   real64 const & oldHydraulicAperture,
                                   real64 const & newHydraulicAperture,
                                   real64 const & dHydraulicAperture_dNormalJump ) const override final
  {
    GEOS_UNUSED_VAR( q );

    compute( oldHydraulicAperture,
             newHydraulicAperture,
             dHydraulicAperture_dNormalJump,
             m_permeability[k][0],
             m_dPerm_dDispJump[k][0] );
  }

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
    GEOS_UNUSED_VAR( dispJump, traction, pressure );

    updateFromAperture( k, q, oldHydraulicAperture, newHydraulicAperture, dHydraulicAperture_dNormalJump );
  }

private:

  arrayView4d< real64 > m_dPerm_dDispJump;
  int m_numDimensionsToUpdate;
};


class ParallelPlatesPermeability : public PermeabilityBase
{
public:

  ParallelPlatesPermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ParallelPlatesPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void initializeState() const override final;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ParallelPlatesPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dDispJump,
                          m_updateTransversalComponent );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * transversalPermeabilityString() { return "transversalPermeability"; }
  } viewKeys;

private:

  array3d< real64 > m_dPerm_dAperture;

  /// Derivative of fracture permeability w.r.t. displacement jump
  array4d< real64 > m_dPerm_dDispJump;

  real64 m_transversalPermeability;

  bool m_updateTransversalComponent;

};

}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_
