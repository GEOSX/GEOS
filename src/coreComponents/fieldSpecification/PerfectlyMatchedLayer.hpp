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

/*
 * PerfectlyMatchedLayer.hpp
 *
 */

#ifndef GEOS_FIELDSPECIFICATION_PERFECTLYMATCHEDLAYER_HPP_
#define GEOS_FIELDSPECIFICATION_PERFECTLYMATCHEDLAYER_HPP_

#include "FieldSpecificationBase.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{

class DomainPartition;

/**
 * @class PerfectlyMatchedLayer
 * A class to manage Perfectly Matched Layer for wave propagation solvers
 */
class PerfectlyMatchedLayer : public FieldSpecificationBase
{
public:
  /**
   * @brief constructor
   * @param name the name of the FieldSpecificationBase in the data repository
   * @param parent the parent group of this group.
   */
  PerfectlyMatchedLayer( string const & name, dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  PerfectlyMatchedLayer() = delete;

  /**
   * @brief destructor
   */
  virtual ~PerfectlyMatchedLayer() = default;

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "PML"; }

  /**
   * @brief Getter for the PML minimum coordinates
   * @return Mininum (x,y,z) coordinates of inner PML boundaries
   */
  R1Tensor32 getMin() const { return m_xMin; }

  /**
   * @brief Getter for the PML maximum coordinates
   * @return Maximum (x,y,z) coordinates of inner PML boundaries
   */
  R1Tensor32 getMax() const { return m_xMax; }

  /**
   * @brief Getter for the PML reflectivity
   * @return Desired reflectivity of the PML region
   */
  real64 getReflectivity() const { return m_reflectivity; }

  /**
   * @brief Getter for the PML thickness
   * @return Thickness of the PML region at left, front, and top sides
   */
  R1Tensor32 getThicknessMinXYZ() const { return m_thicknessMinXYZ; }

  /**
   * @brief Getter for the PML thickness
   * @return Thickness of the PML region at right, back, and bottom sides
   */
  R1Tensor32 getThicknessMaxXYZ() const { return m_thicknessMaxXYZ; }

  /**
   * @brief Getter for the PML wave speed
   * @return Wave speed of the PML region at left, front, and top sides
   */
  R1Tensor32 getWaveSpeedMinXYZ() const { return m_waveSpeedMinXYZ; }

  /**
   * @brief Getter for the PML wave speed
   * @return Wave speed of the PML region at right, back, and bottom sides
   */
  R1Tensor32 getWaveSpeedMaxXYZ() const { return m_waveSpeedMaxXYZ; }

  /**
   * @brief View keys
   */
  struct viewKeyStruct : public FieldSpecificationBase::viewKeyStruct
  {
    /// @return String key for the mininum (x,y,z) coordinates of inner PML boundaries
    static constexpr char const * xMinString() { return "xMin"; }

    /// @return String key for the maximum (x,y,z) coordinates of inner PML boundaries
    static constexpr char const * xMaxString() { return "xMax"; }

    /// @return String key for the reflectivity of the PML region
    static constexpr char const * reflectivityString() { return "reflectivity"; }

    /// @return String key for the thickness of the PML region (left, front, top sides)
    static constexpr char const * thicknessMinXYZString() { return "thicknessMinXYZ"; }

    /// @return String key for the thickness of the PML region (right, back, bottom sides)
    static constexpr char const * thicknessMaxXYZString() { return "thicknessMaxXYZ"; }

    /// @return String key for the wave speed in the PML region (left, front, top sides)
    static constexpr char const * waveSpeedMinXYZString() { return "waveSpeedMinXYZ"; }

    /// @return String key for the wave speed in the PML region (right, back, bottom sides)
    static constexpr char const * waveSpeedMaxXYZString() { return "waveSpeedMaxXYZ"; }

  };

  /**
   * @brief Safeguard for the minimum allowed PML thickness
   */
  static constexpr real64 minThickness = 0.001;

  /**
   * @brief Smallest possible values for xMin, below which they are computed internally
   */
  static constexpr real64 smallestXMin = -999999.0;

  /**
   * @brief Largest possible values for xMax, below which they are computed internally
   */
  static constexpr real64 largestXMax = 999999.0;

protected:

  virtual void postInputInitialization() override final;

private:

  /// Mininum (x,y,z) coordinates of inner PML boundaries
  R1Tensor32 m_xMin;

  /// Maximum (x,y,z) coordinates of inner PML boundaries
  R1Tensor32 m_xMax;

  /// Desired reflectivity of the PML region, used to compute the damping profile
  real32 m_reflectivity;

  /// Thickness of the PML region, used to compute the damping profile
  R1Tensor32 m_thicknessMinXYZ;
  R1Tensor32 m_thicknessMaxXYZ;

  /// Wave speed in the PML region, used to compute the damping profile
  R1Tensor32 m_waveSpeedMinXYZ;
  R1Tensor32 m_waveSpeedMaxXYZ;

};

} /* namespace geos */

#endif /* GEOS_FIELDSPECIFICATION_PERFECTLYMATCHEDLAYER_HPP_ */
