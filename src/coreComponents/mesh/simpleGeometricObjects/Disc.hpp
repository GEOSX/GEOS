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
 * @file Disc.hpp
 */

#ifndef GEOS_MESH_SIMPLEGEOMETRICOBJECTS_DISC_HPP_
#define GEOS_MESH_SIMPLEGEOMETRICOBJECTS_DISC_HPP_

#include "PlanarGeometricObject.hpp"

namespace geos
{

/**
 * @class Disc
 * @brief Class to represent a geometric disc in GEOSX.
 */
class Disc : public PlanarGeometricObject
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  Disc( const string & name,
        Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~Disc() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "Disc"; }

  ///@}

  bool isCoordInObject( real64 const ( &coord ) [3] ) const override final;

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get the center of the disc.
   * @return the center of the disc
   */
  virtual R1Tensor & getCenter() override final {return m_center;}

  /**
   * @copydoc getCenter()
   */
  virtual R1Tensor const & getCenter() const override final {return m_center;}

protected:

  /**
   * @brief This function provides the capability to post process input values prior to
   * any other initialization operations.
   */
  virtual void postInputInitialization() override final;

private:

  /// center of the disc in (x,y,z) coordinates
  R1Tensor m_center;
  /// Dimensions of the disc's radius
  real64 m_radius;
  /// tolerance to determine if a point sits on the disc or not
  real64 m_tolerance;

  /// @cond DO_NOT_DOCUMENT

  struct viewKeyStruct
  {
    static constexpr char const * centerString() { return "center"; }
    static constexpr char const * radiusString() { return "radius"; }
    static constexpr char const * toleranceString() { return "tolerance"; }
  };

  /// @endcond

};
} /* namespace geos */

#endif /* GEOS_MESH_SIMPLEGEOMETRICOBJECTS_DISC_HPP_*/
