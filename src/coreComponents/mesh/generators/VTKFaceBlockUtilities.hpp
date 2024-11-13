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

#ifndef GEOS_VTKFACEBLOCKUTILITIES_HPP
#define GEOS_VTKFACEBLOCKUTILITIES_HPP

#include "CellBlockManager.hpp"

#include "common/DataTypes.hpp"

#include <vtkDataSet.h>
#include <vtkSmartPointer.h>

namespace geos::vtk
{

/**
 * @brief Attach the face block information to the cell block manager.
 * @param faceBlockName[in] The name of the face block.
 * @param faceMesh[in] The vtk mesh for the face block.
 * @param mesh[in] The vtk volumic mesh.
 * @param cellBlockManager[inout] The cell block manager that will receive the face block information.
 */
void importFractureNetwork( string const & faceBlockName,
                            vtkSmartPointer< vtkDataSet > faceMesh,
                            vtkSmartPointer< vtkDataSet > mesh,
                            CellBlockManager & cellBlockManager );
}

#endif // include guard
