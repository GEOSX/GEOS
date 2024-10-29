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
 * @file VTKWellGenerator.hpp
 */

#ifndef GEOS_MESH_GENERATORS_VTKWELLGENERATOR_HPP
#define GEOS_MESH_GENERATORS_VTKWELLGENERATOR_HPP

#include "WellGeneratorBase.hpp"

#include "VTKUtilities.hpp"

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"

#include <vtkDataSet.h>


namespace geos
{

/**
 *  @class VTKWellGenerator
 *  @brief The VTKWellGenerator class provides a class implementation of VTK generated well.
 */
class VTKWellGenerator final : public WellGeneratorBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Main constructor for VTKWellGenerator base class.
   * @param[in] name of the VTKWellGenerator object
   * @param[in] parent the parent Group pointer for the WellGenerator object
   */
  VTKWellGenerator( const string & name,
                    Group * const parent );

/**
 * @brief Return the name of the VTKWellGenerator in object Catalog.
 * @return string that contains the key name to VTKWellGenerator in the Catalog
 */
  static string catalogName() { return "VTKWell"; }

  ///@}


private:

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * filePathString() { return "file"; }
  };
  /// @endcond


  /**
   * @brief
   */
  void fillPolylineDataStructure( ) override;


  /// Path to the mesh file
  Path m_filePath;

};

} // namespace geos

#endif /* GEOS_MESH_GENERATORS_VTKWELLGENERATOR_HPP */
