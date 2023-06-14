/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "RelpermDriverRunTest.hpp"
#include "BrooksCoreyStone2RelativePermeability.hpp"


namespace geos
{
    template void RelpermDriver::runTest< geosx::constitutive::BrooksCoreyStone2RelativePermeability >( geosx::constitutive::BrooksCoreyStone2RelativePermeability &, arrayView2d< real64 > const & );
}
