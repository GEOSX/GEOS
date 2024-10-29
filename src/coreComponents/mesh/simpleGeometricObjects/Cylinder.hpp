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
 * @file Cylinder.hpp
 *
 */

#ifndef GEOS_MESH_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_
#define GEOS_MESH_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geos
{

/**
 * @class Cylinder
 * @brief Class to represent a geometric cylinder in GEOSX.
 */
class Cylinder : public SimpleGeometricObjectBase
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
  Cylinder( const string & name,
            Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~Cylinder() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "Cylinder"; }

  ///@}

  bool isCoordInObject( real64 const ( &coord ) [3] ) const override final;

  /// @cond DO_NOT_DOCUMENT

  struct viewKeyStruct
  {
    static constexpr char const * point1String() { return "firstFaceCenter"; }
    static constexpr char const * point2String() { return "secondFaceCenter"; }
    static constexpr char const * radiusString() { return "outerRadius"; }
    static constexpr char const * innerRadiusString() { return "innerRadius"; }
  };

  /// @endcond

private:

  /// Center point of the first face of the cylinder
  R1Tensor m_point1;

  /// Center point of the second face of the cylinder
  R1Tensor m_point2;

  /// Outer radius of the cylinder
  real64 m_radius = 0.0;

  /// Inner radius of the annulus
  real64 m_innerRadius = 0.0;

};
} /* namespace geos */

#endif /* GEOS_MESH_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_
        */
