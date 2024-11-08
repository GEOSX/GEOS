/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_VTKEMBEDDEDSURFACEBLOCKUTILITIES_HPP
#define GEOS_VTKEMBEDDEDSURFACEBLOCKUTILITIES_HPP

#include "CellBlockManager.hpp"

#include "common/DataTypes.hpp"

#include <vtkDataSet.h>
#include <vtkSmartPointer.h>

namespace geos::vtk
{

/**
 * @brief Attach the embedded surface block information to the cell block manager.
 * @param embeddedSurfaceBlockName[in] The name of the embedded surface block.
 * @param embeddedSurfaceMesh[in] The vtk mesh for the embedded surface block.
 * @param mesh[in] The vtk volumic mesh.
 * @param cellBlockManager[inout] The cell block manager that will receive the embedded surface block information.
 */
void importEmbeddedFractureNetwork( string const & embeddedSurfaceBlockName,
                                    vtkSmartPointer< vtkDataSet > embeddedSurfaceMesh,
                                    vtkSmartPointer< vtkDataSet > mesh,
                                    CellBlockManager & cellBlockManager );
}

#endif // include guard
