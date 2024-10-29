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
 * @file TableRelativePermeabilityHelpers.hpp
 */

#ifndef GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_TABLERELATIVEPERMEABILITYHELPERS_HPP
#define GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_TABLERELATIVEPERMEABILITYHELPERS_HPP

#include "functions/TableFunction.hpp"

namespace geos
{
namespace constitutive
{

struct TableRelativePermeabilityHelpers
{

  /**
   * @brief Validate the relative permeability table provided in input (increasing phase vol frac and rel perm, etc)
   * @param[in] relPermTable the relative permeability table (kr vs s) for a given phase)
   * @param[in] fullConstitutiveName name of the constitutive model
   * @param[out] phaseMinVolFrac the phase minimum volume fraction read from the table
   * @param[out] phaseMaxVolFrac the phase maximum volume fraction read from the table
   * @param[out] phaseRelPermEndPoint the end-point relative permeability
   */
  static
  void validateRelativePermeabilityTable( TableFunction const & relPermTable,
                                          string const & fullConstitutiveName,
                                          real64 & phaseMinVolFrac,
                                          real64 & phaseMaxVolFrac,
                                          real64 & phaseRelPermEndPoint );

};

} // namespace constitutive

} // namespace geos

#endif // GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_TABLECAPILLARYPRESSUREHELPERS_HPP
