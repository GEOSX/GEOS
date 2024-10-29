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

/*
 * @file Box.hpp
 *
 */

#ifndef GEOS_MESH_SIMPLEGEOMETRICOBJECTS_BOX_HPP_
#define GEOS_MESH_SIMPLEGEOMETRICOBJECTS_BOX_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geos
{

/**
 * @class Box
 * @brief Class to represent a geometric box in GEOSX.
 */
class Box : public SimpleGeometricObjectBase
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
  Box( const string & name,
       Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~Box() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "Box"; }

  ///@}

  bool isCoordInObject( real64 const ( &coord ) [3] ) const override final;

protected:

  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  virtual void postInputInitialization() override final;

private:

  /// Mininum (x,y,z) coordinates of the box
  R1Tensor m_min;
  /// Maximum (x,y,z) coordinates of the box
  R1Tensor m_max;
  /// Strike angle of the box
  real64 m_strikeAngle=0.0;
  /// Coordinates of the center of the box
  R1Tensor m_boxCenter={0.0, 0.0, 0.0};
  /// Cosine of the strike angle of the box
  real64 m_cosStrike=0.0;
  /// Sine of the strike angle of the box
  real64 m_sinStrike=0.0;

  /// @cond DO_NOT_DOCUMENT

  struct viewKeyStruct
  {
    static constexpr char const * xMinString() { return "xMin"; }
    static constexpr char const * xMaxString() { return "xMax"; }
    static constexpr char const * strikeAngleString() { return "strike"; }
    static constexpr char const * boxCenterString() { return "center"; }
    static constexpr char const * cosStrikeString() { return "cosStrike"; }
    static constexpr char const * sinStrikeString() { return "sinStrike"; }

    dataRepository::ViewKey xmin = { xMinString() };
    dataRepository::ViewKey xmax = { xMaxString() };
  } viewKeys;

  /// @endcond

};
} /* namespace geos */

#endif /* GEOS_MESH_SIMPLEGEOMETRICOBJECTS_BOX_HPP_ */
