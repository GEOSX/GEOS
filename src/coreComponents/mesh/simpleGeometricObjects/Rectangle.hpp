/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Rectangle.hpp
 */

#ifndef GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_RECTANGLE_HPP_
#define GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_RECTANGLE_HPP_

#include "PlanarGeometricObject.hpp"

namespace geos
{

/**
 * @class Rectangle
 * @brief Class to represent a geometric box in GEOSX.
 */
class Rectangle : public PlanarGeometricObject
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
  Rectangle( const string & name,
             Group * const parent );

  /**
   * @brief Internal constructor. This is used to make planar cuts from point (oldX, oldY) to (newX, newY)
   * in 2.5D problems.
   * @param oldX x-coordinate of first point
   * @param oldY y-coordinate of first point
   * @param newX x-coordinate of second point
   * @param newY y-coordinate of second point
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  Rectangle( const real64 oldX, const real64 oldY, const real64 newX,
             const real64 newY, const string & name, Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~Rectangle() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "Rectangle"; }

  ///@}

  bool isCoordInObject( real64 const ( &coord ) [3] ) const override final;

  /**
   * @brief Find the bounds of the plane.
   */
  void findRectangleLimits();

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get the origin of the plane.
   * @return the origin of the plane
   */
  virtual R1Tensor & getCenter() override final {return m_origin;}

  /**
   * @copydoc getCenter()
   */
  virtual R1Tensor const & getCenter() const override final {return m_origin;}



protected:

  /**
   * @brief This function provides the capability to post process input values prior to
   * any other initialization operations.
   */
  virtual void postInputInitialization() override final;

private:

  /// Origin point (x,y,z) of the plane (basically, any point on the plane)
  R1Tensor m_origin;
  /// Dimensions of the bounded plane
  array1d< real64 > m_dimensions;
  /// Length and width of the bounded plane
  array2d< real64 > m_points;
  /// tolerance to determine if a point sits on the plane or not
  real64 m_tolerance;

  /// @cond DO_NOT_DOCUMENT

  struct viewKeyStruct
  {
    static constexpr char const * originString() { return "origin"; }
    static constexpr char const * dimensionsString() { return "dimensions"; }
    static constexpr char const * toleranceString() { return "tolerance"; }
  };

  /// @endcond

};
} /* namespace geosx */

#endif /* GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_RECTANGLE_HPP_*/
