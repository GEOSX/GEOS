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
 * PVTDriverRunTestCompositionalMultiphaseFluid.cpp
 */

#include "constitutive/fluid/multifluid/PVTDriverRunTest.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluidPVTPackage.hpp"

namespace geos
{
template void PVTDriver::runTest< constitutive::CompositionalMultiphaseFluidPVTPackage >( constitutive::CompositionalMultiphaseFluidPVTPackage &, arrayView2d< real64 > const & );
}
