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
 * @file LinearIsotropicDispersion.hpp
 */

#ifndef GEOS_CONSTITUTIVE_DISPERSION_LINEARISOTROPICDISPERSION_HPP_
#define GEOS_CONSTITUTIVE_DISPERSION_LINEARISOTROPICDISPERSION_HPP_

#include "constitutive/dispersion/DispersionBase.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief The update class for linear isotropic diffusion
 */
class LinearIsotropicDispersionUpdate : public DispersionBaseUpdate
{
public:

  /**
   * @brief Constructor for the class performing the dispersion updates
   * @param dispersivity the array of cell-wise dispersivities in the subregion
   * @param longitunidalDispersivity longitudinal dispersivity in the subregion
   */
  LinearIsotropicDispersionUpdate( arrayView3d< real64 > const & dispersivity,
                                   real64 const & longitudinalDispersivity )
    : DispersionBaseUpdate( dispersivity ),
    m_longitudinalDispersivity( longitudinalDispersivity )
  {}

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice1d< real64 const > const & laggedTotalVelocityComponents ) const override
  {
    real64 const velocityNorm =
      sqrt( laggedTotalVelocityComponents[0]*laggedTotalVelocityComponents[0]
            + laggedTotalVelocityComponents[1]*laggedTotalVelocityComponents[1]
            + laggedTotalVelocityComponents[2]*laggedTotalVelocityComponents[2] );
    for( integer i = 0; i < 3; ++i )
    {
      m_dispersivity[k][q][i] = m_longitudinalDispersivity * velocityNorm;
    }
  }

protected:

  /// Longitudinal dispersivity
  real64 const m_longitudinalDispersivity;

};

/**
 * @brief The class for constant dispersion
 */
class LinearIsotropicDispersion : public DispersionBase
{
public:

  /**
   * @brief Constructor for the class storing constant dispersion
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  LinearIsotropicDispersion( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  static string catalogName() { return "LinearIsotropicDispersion"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void initializeVelocityState( arrayView2d< real64 const > const & initialVelocity ) const override;

  virtual void saveConvergedVelocityState( arrayView2d< real64 const > const & convergedVelocity ) const override;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = LinearIsotropicDispersionUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_dispersivity,
                          m_longitudinalDispersivity );
  }

  struct viewKeyStruct : public DispersionBase::viewKeyStruct
  {
    static constexpr char const * longitudinalDispersivityString() { return "longitudinalDispersivity"; }
  } viewKeys;

protected:

  virtual void postInputInitialization() override;

private:

  /// Longitudinal dispersivity
  real64 m_longitudinalDispersivity;

};

} // namespace constitutive

} // namespace geos


#endif // GEOS_CONSTITUTIVE_DISPERSION_LINEARISOTROPICDISPERSION_HPP_
