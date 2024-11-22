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
 * @file WaveSolverTypeDefDG.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERTYPEDEFDG_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERTYPEDEFDG_HPP_


#if !defined( GEOS_USE_HIP )
#include "finiteElement/elementFormulations/BB_Tetrahedron.hpp"
#endif


#if !defined( GEOS_USE_HIP )
#define DG_FE_TYPES \
  finiteElement::BB1_Tetrahedron, \
  finiteElement::BB2_Tetrahedron, \
  finiteElement::BB3_Tetrahedron, \
  finiteElement::BB4_Tetrahedron, \
  finiteElement::BB5_Tetrahedron
#else                                                                           
#define DG_FE_TYPES                                                            
#endif                                                                          


#define SELECTED_FE_TYPES DG_FE_TYPES

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERTYPEDEFDG_HPP_ */
