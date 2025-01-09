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

/***
 *  @file KilloughHysteresis.hpp
 */

#include "KilloughHysteresis.hpp"


namespace geos
{

using namespace dataRepository;

namespace constitutive
{


void KilloughHysteresis::postProcessInput( real64 const & jerauldParam_a,
                                           real64 const & jerauldParam_b,
                                           real64 const & killoughCurvatureParamRelPerm )
{
  GEOS_THROW_IF( jerauldParam_a < 0,
                 GEOS_FMT( "{}: the parameter {} must be positive",
                           catalogName(),
                           viewKeyStruct::jerauldParameterAString() ),
                 InputError );

  GEOS_THROW_IF( jerauldParam_b < 0,
                 GEOS_FMT( "{}: the paramater {} must be postitive",
                           catalogName(),
                           viewKeyStruct::jerauldParameterBString() ),
                 InputError );

  GEOS_THROW_IF( killoughCurvatureParamRelPerm < 0,
                 GEOS_FMT( "{}: the paramater {} must be postitive",
                           catalogName(),
                           viewKeyStruct::killoughCurvatureParameterString() ),
                 InputError );

}



} //end namespace

} //end namespace
