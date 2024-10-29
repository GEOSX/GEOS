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
 * @file PlanarGeometricObject.hpp
 */

#ifndef GEOS_MESH_SIMPLEGEOMETRICOBJECTS_PLANARGEOMETRICOBJECT_HPP_
#define GEOS_MESH_SIMPLEGEOMETRICOBJECTS_PLANARGEOMETRICOBJECT_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geos
{

/**
 * @class PlanarGeometricObject
 * @brief Abstract class to implement functions used by all bounded geometric objects in GEOSX, such as disc or plane.
 */
class PlanarGeometricObject : public SimpleGeometricObjectBase
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
  PlanarGeometricObject( const string & name,
                         Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~PlanarGeometricObject() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "PlanarGeometricObject"; }

  ///@}

  /**
   * @brief Check if the input coordinates are in the object.
   * @param[in] coord the coordinates to test
   * @return true if the coordinates are in the object, false otherwise
   */
  virtual bool isCoordInObject( real64 const ( &coord ) [3] ) const override = 0;

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get the normal to the plane.
   * @return the normal vector
   */
  R1Tensor & getNormal() {return m_normal;}

  /**
   * @copydoc getNormal()
   */
  R1Tensor const & getNormal() const {return m_normal;}

  /**
   * @brief Get one of the tangent vectors defining the orthonormal basis along with the normal.
   * @return the tangent vector
   */
  R1Tensor & getWidthVector() {return m_widthVector;}

  /**
   * @copydoc getWidthVector()
   */
  R1Tensor const & getWidthVector() const {return m_widthVector;}

  /**
   * @brief Get one of the tangent vectors defining the orthonormal basis along with the normal.
   * @return the length vector
   */
  R1Tensor & getLengthVector() {return m_lengthVector;}

  /**
   * @copydoc getLengthVector()
   */
  R1Tensor const & getLengthVector() const {return m_lengthVector;}

  /**
   * @brief Get the origin of the plane.
   * @return the origin of the plane
   */
  virtual R1Tensor & getCenter() = 0;

  /**
   * @copydoc getCenter()
   */
  virtual R1Tensor const & getCenter() const = 0;

protected:

  /// Normal (n_x,n_y,n_z) to the plane (will be normalized automatically)
  R1Tensor m_normal;
  /// Length vector in the orthonormal basis along with the normal
  R1Tensor m_lengthVector;
  /// Width vector in the orthonormal basis along with the normal
  R1Tensor m_widthVector;
  /// tolerance to determine if a point sits on the PlanarGeometricObject or not
  real64 m_tolerance;

  /// tolerance to check if base is orthonormal
  static constexpr real64 orthoNormalBaseTolerance = 1e-10;

  /// @cond DO_NOT_DOCUMENT

  struct viewKeyStruct
  {
    static constexpr char const * normalString() { return "normal"; }
    static constexpr char const * mLengthVectorString() { return "lengthVector"; }
    static constexpr char const * mWidthVectorString() { return "widthVector"; }
  };

  /// @endcond

};
} /* namespace geos */

#endif /* GEOS_MESH_SIMPLEGEOMETRICOBJECTS_PLANARGEOMETRICOBJECT_HPP_*/
