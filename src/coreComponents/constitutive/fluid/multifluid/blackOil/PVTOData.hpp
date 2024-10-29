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
 * @file PVTOData.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_BLACKOIL_PVTODATA_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_BLACKOIL_PVTODATA_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

struct PVTOData
{
public:

  class KernelWrapper
  {
public:

    /// @cond DO_NOT_DOCUMENT
    /// We need these SMFs to avoid host-device errors with CUDA.
    KernelWrapper() = default;
    KernelWrapper( KernelWrapper const & ) = default;
    KernelWrapper & operator=( KernelWrapper const & ) = default;
    KernelWrapper & operator=( KernelWrapper && ) = default;
    /// @endcond

    /**
     * @brief The constructor of the kernel wrapper
     */
    KernelWrapper( arrayView1d< real64 const > const & Rs,
                   arrayView1d< real64 const > const & bubblePressure,
                   arrayView1d< real64 const > const & saturatedBo,
                   arrayView1d< real64 const > const & saturatedViscosity,
                   arrayView2d< real64 const > const & undersaturatedPressure,
                   arrayView2d< real64 const > const & undersaturatedBo,
                   arrayView2d< real64 const > const & undersaturatedViscosity,
                   arrayView1d< real64 const > const & surfaceMassDensity,
                   arrayView1d< real64 const > const & surfaceMoleDensity );

    void move( LvArray::MemorySpace const space, bool const touch ) const
    {
      m_Rs.move( space, touch );
      m_bubblePressure.move( space, touch );

      m_saturatedBo.move( space, touch );
      m_saturatedViscosity.move( space, touch );

      m_undersaturatedPressure2d.move( space, touch );
      m_undersaturatedBo2d.move( space, touch );
      m_undersaturatedViscosity2d.move( space, touch );

      m_surfaceMassDensity.move( space, touch );
      m_surfaceMoleDensity.move( space, touch );
    }

public:

    /// Rs array
    arrayView1d< real64 const > m_Rs;
    /// Bubble-point pressure array
    arrayView1d< real64 const > m_bubblePressure;

    // Saturated data (free gas phase present)

    /// Saturated oil phase formation volume factor
    arrayView1d< real64 const > m_saturatedBo;
    /// Saturated oil phase viscosity
    arrayView1d< real64 const > m_saturatedViscosity;

    // Undersaturated data (no free gas)

    /// Undersaturated pressure
    arrayView2d< real64 const > m_undersaturatedPressure2d;
    /// Undersaturated oil phase formation volume factor
    arrayView2d< real64 const > m_undersaturatedBo2d;
    /// Undersaturated oil phase viscosity
    arrayView2d< real64 const > m_undersaturatedViscosity2d;

    /// Surface mass density
    arrayView1d< real64 const > m_surfaceMassDensity;
    /// Surface mole density
    arrayView1d< real64 const > m_surfaceMoleDensity;

  };

public:

  KernelWrapper createKernelWrapper() const;

  /// Max pressure used for the construction of the undersaturated tables
  real64 maxRelativePressure = 0.0;

  /// Rs array
  array1d< real64 > Rs;
  /// Bubble-point pressure array
  array1d< real64 > bubblePressure;

  // Saturated data (free gas phase present)

  /// Number of saturated points
  integer numSaturatedPoints = 0;
  /// Saturated oil phase formation volume factor
  array1d< real64 > saturatedBo;
  /// Saturated oil phase viscosity
  array1d< real64 > saturatedViscosity;

  // Undersaturated data (no free gas)

  /// Undersaturated pressure
  array2d< real64 > undersaturatedPressure2d;
  /// Undersaturated oil phase formation volume factor
  array2d< real64 > undersaturatedBo2d;
  /// Undersaturated oil phase viscosity
  array2d< real64 > undersaturatedViscosity2d;

  /// Surface mass density
  array1d< real64 > surfaceMassDensity;
  /// Surface mole density
  array1d< real64 > surfaceMoleDensity;

  // temporary arrays used in the construction of the undersaturated arrays
  // can be discarded after construction
  array1d< array1d< real64 > > undersaturatedPressure;
  array1d< array1d< real64 > > undersaturatedBo;
  array1d< array1d< real64 > > undersaturatedViscosity;

};

} //namespace constitutive

} //namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTODATA_HPP_
