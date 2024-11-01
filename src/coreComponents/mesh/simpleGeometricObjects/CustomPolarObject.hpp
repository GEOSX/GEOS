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
 * @file CustomPolarObject.hpp
 */

#ifndef GEOS_MESH_SIMPLEGEOMETRICOBJECTS_CUSTOMPOLAROBJECT_HPP_
#define GEOS_MESH_SIMPLEGEOMETRICOBJECTS_CUSTOMPOLAROBJECT_HPP_

#include "PlanarGeometricObject.hpp"

namespace geos
{

/**
 * @class CustomPolarObject
 * @brief Class to represent a geometric disc in GEOSX.
 */
class CustomPolarObject : public PlanarGeometricObject
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
  CustomPolarObject( const string & name,
                     Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~CustomPolarObject() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "CustomPolarObject"; }

  ///@}

  bool isCoordInObject( real64 const ( &coord ) [3] ) const override final;

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get the center of the CustomPolarObject.
   * @return the center of the CustomPolarObject
   */
  virtual R1Tensor & getCenter() override final {return m_center;}

  /**
   * @copydoc getCenter()
   */
  virtual R1Tensor const & getCenter() const override final {return m_center;}

  /**
   * @brief Update the geometric function describing the boundary of the object.
   * @param coefficients define all the coefficients of the radius function.
   */
  void setCustomPolarObjectFunction( array1d< real64 > & coefficients )
  {
    //m_radius.m_coefficients = coefficients;
    m_coefficients = coefficients;
  }

  /**
   * @brief Get value of the radius of the surface for each angle theta
   * @param angle the given angle
   * @return the radius for that angle
   */
  real64 getRadius( real64 angle ) const
  {
    real64 radius = 0;
    integer count = 0;
    for( auto coeff:m_coefficients )
    {
      radius = radius + coeff*cos( count*angle );
      count++;
    }
    return radius;
  }

protected:

  /**
   * @brief This function provides the capability to post process input values prior to
   * any other initialization operations.
   */
  virtual void postInputInitialization() override final;

private:

  /// center of the CustomPolarObject in (x,y,z) coordinates
  R1Tensor m_center;
  /// Coefficients of the polar function relating the radius to theta
  array1d< real64 > m_coefficients;
  /// tolerance to determine if a point sits on the CustomPolarObject or not
  real64 m_tolerance;

  /// @cond DO_NOT_DOCUMENT

  struct viewKeyStruct
  {
    static constexpr char const * centerString() { return "center"; }
    static constexpr char const * coefficientsString() { return "coefficients"; }
    static constexpr char const * toleranceString() { return "tolerance"; }
  };

  /// @endcond

};
} /* namespace geos */

#endif /* GEOS_MESH_SIMPLEGEOMETRICOBJECTS_CustomPolarObject_HPP_*/
