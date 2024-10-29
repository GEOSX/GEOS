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
 * @file TableRelativePermeabilityHelpers.cpp
 */

#include "TableRelativePermeabilityHelpers.hpp"

#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

void
TableRelativePermeabilityHelpers::validateRelativePermeabilityTable( TableFunction const & relPermTable,
                                                                     string const & fullConstitutiveName,
                                                                     real64 & phaseMinVolFrac,
                                                                     real64 & phaseMaxVolFrac,
                                                                     real64 & phaseRelPermEndPoint )
{
  ArrayOfArraysView< real64 const > coords = relPermTable.getCoordinates();

  GEOS_THROW_IF_NE_MSG( relPermTable.getInterpolationMethod(), TableFunction::InterpolationType::Linear,
                        GEOS_FMT( "{}: TableFunction '{}' interpolation method must be linear",
                                  fullConstitutiveName, relPermTable.getDataContext() ),
                        InputError );
  GEOS_THROW_IF_NE_MSG( relPermTable.numDimensions(), 1,
                        GEOS_FMT( "{}: TableFunction '{}' must have a single independent coordinate",
                                  fullConstitutiveName, relPermTable.getDataContext() ),
                        InputError );
  GEOS_THROW_IF_LT_MSG( coords.sizeOfArray( 0 ), 2,
                        GEOS_FMT( "{}: TableFunction `{}` must contain at least two values",
                                  fullConstitutiveName, relPermTable.getDataContext() ),
                        InputError );

  arraySlice1d< real64 const > phaseVolFrac = coords[0];
  arrayView1d< real64 const > const relPerm = relPermTable.getValues();
  phaseMinVolFrac = phaseVolFrac[0];
  phaseMaxVolFrac = phaseVolFrac[phaseVolFrac.size()-1];
  phaseRelPermEndPoint = relPerm[relPerm.size()-1];

  // note that the TableFunction class has already checked that coords.sizeOfArray( 0 ) == relPerm.size()
  GEOS_THROW_IF( !isZero( relPerm[0] ),
                 GEOS_FMT( "{}: TableFunction '{}' first value must be equal to 0",
                           fullConstitutiveName, relPermTable.getDataContext() ),
                 InputError );
  for( localIndex i = 1; i < coords.sizeOfArray( 0 ); ++i )
  {
    // check phase volume fraction
    GEOS_THROW_IF( phaseVolFrac[i] < 0 || phaseVolFrac[i] > 1,
                   GEOS_FMT( "{}: TableFunction '{}' values must be between 0 and 1",
                             fullConstitutiveName, relPermTable.getDataContext() ),
                   InputError );

    // note that the TableFunction class has already checked that the coordinates are monotone

    // check phase relative permeability
    GEOS_THROW_IF( !isZero( relPerm[i] ) && (relPerm[i] - relPerm[i-1]) < 1e-15,
                   GEOS_FMT( "{}: TableFunction '{}' values must be strictly increasing (|Delta kr| > 1e-15 between two non-zero values)",
                             fullConstitutiveName, relPermTable.getDataContext() ),
                   InputError );

    if( isZero( relPerm[i-1] ) && !isZero( relPerm[i] ) )
    {
      phaseMinVolFrac = phaseVolFrac[i-1];
    }
  }
}


} // namespace constitutive

} // namespace geos
