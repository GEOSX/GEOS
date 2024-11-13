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
 * @file DispersionBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_DISPERSION_DISPERSIONBASE_HPP
#define GEOS_CONSTITUTIVE_DISPERSION_DISPERSIONBASE_HPP

#include "common/DataLayouts.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief The abstract base class to perform the dispersion update
 */
class DispersionBaseUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_dispersivity.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_dispersivity.size( 1 ); }

protected:

  /**
   * @brief Constructor for the class performing the dispersion updates
   * @param dispersivity the array of cell-wise dispersion in the subregion
   */
  DispersionBaseUpdate( arrayView3d< real64 > const & dispersivity )
    : m_dispersivity( dispersivity )
  {}

  /// View on the cell-wise dispersivity
  arrayView3d< real64 > const m_dispersivity;

private:

  /**
   * @brief Pointwise update function called from the solver
   * @param[in] k index of the cell in the subRegion
   * @param[in] q constitutive index (equal to one in this class)
   * @param[in] laggedTotalVelocityComponents the components of the total velocity
   */
  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const > const & laggedTotalVelocityComponents ) const = 0;
};

/**
 * @brief The abstract base class for dispersion
 */
class DispersionBase : public ConstitutiveBase
{
public:

  /**
   * @brief Constructor for the abstract base class
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  DispersionBase( string const & name, dataRepository::Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @brief Getter for the dispersivities in the subRegion
   * @return an arrayView of dispersivities
   */
  arrayView3d< real64 const > dispersivity() const { return m_dispersivity; }

  /**
   * @brief Initialize the velocity state (needed because dispersion depends on total velocity)
   * @param[in] initialVelocity the initial velocity field after reservoir initialization
   *
   * Note: this is needed because for now, the velocity field is treated **explicitly** in the dispersion tensor
   */
  virtual void initializeVelocityState( arrayView2d< real64 const > const & initialVelocity ) const
  { GEOS_UNUSED_VAR( initialVelocity ); }

  /**
   * @brief Save the velocity state (needed because dispersion depends on total velocity)
   * @param[in] convergedVelocity the converged velocity field
   *
   * Note: this is needed because for now, the velocity is treated **explicitly** in the dispersion tensor
   */
  virtual void saveConvergedVelocityState( arrayView2d< real64 const > const & convergedVelocity ) const
  { GEOS_UNUSED_VAR( convergedVelocity ); }

private:

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  void resizeFields( localIndex const size, localIndex const numPts );

protected:

  virtual void postInputInitialization() override;

  /// cell-wise dispersivity in the subregion
  /// TODO: support full tensor if linear isotropic diffusion is no longer enough
  array3d< real64 > m_dispersivity;

};

} // namespace constitutive

} // namespace geos


#endif // GEOS_CONSTITUTIVE_DISPERSION_DISPERSIONBASE_HPP
