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
 * @file PerforationFields.hpp
 */

#ifndef GEOS_MESH_PERFORATIONFIELDS_HPP_
#define GEOS_MESH_PERFORATIONFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace perforation
{

///@cond DO_NOT_DOCUMENT
DECLARE_FIELD( reservoirElementRegion,
               "reservoirElementRegion",
               array1d< localIndex >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "For each perforation, elementRegion index of the perforated element" );

DECLARE_FIELD( reservoirElementSubRegion,
               "reservoirElementSubregion",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "For each perforation, elementSubRegion index of the perforated element" );

DECLARE_FIELD( reservoirElementIndex,
               "reservoirElementIndex",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "For each perforation, element index of the perforated element" );

DECLARE_FIELD( wellElementIndex,
               "wellElementIndex",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "For each perforation, index of the well element" );

DECLARE_FIELD( wellTransmissibility,
               "wellTransmissibility",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "For each perforation, well transmissibility" );

DECLARE_FIELD( wellSkinFactor,
               "wellSkinFactor",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "For each perforation, well skin factor" );

DECLARE_FIELD( location,
               "location",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "For each perforation, physical location (x,y,z coordinates)" );

///@endcond DO_NOT_DOCUMENT
}

}

}

#endif // GEOS_MESH_PERFORATIONFIELDS_HPP_
