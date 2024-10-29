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

#include "constitutive/fluid/multifluid/PVTDriverRunTest.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/CO2BrineFluid.hpp"

namespace geos
{
template void PVTDriver::runTest< constitutive::CO2BrinePhillipsFluid >( constitutive::CO2BrinePhillipsFluid &, arrayView2d< real64 > const & );
}
