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

/*
 * SourceFluxBoundaryCondition.cpp
 *
 */

#include "SourceFluxBoundaryCondition.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"

namespace geos
{
using namespace dataRepository;

SourceFluxBoundaryCondition::SourceFluxBoundaryCondition( string const & name, Group * const parent ):
  FieldSpecificationBase( name, parent )
{
  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).
    setInputFlag( InputFlags::FALSE );
  setFieldName( catalogName() );

  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::functionNameString() ).
    setDescription( GEOS_FMT( "Name of a function that specifies the variation of the production rate variations of this {}."
                              "Multiplied by {}. If no function is provided, a constant value of 1 is used."
                              "The producted fluid rate unit is in kg by default, or in mole if the flow solver has a {} of 0.",
                              catalogName(),
                              FieldSpecificationBase::viewKeyStruct::scaleString(),
                              CompositionalMultiphaseBase::viewKeyStruct::useMassFlagString() ) );

  getWrapper< real64 >( FieldSpecificationBase::viewKeyStruct::scaleString() ).
    setDescription( GEOS_FMT( "Multiplier of the {0} value. If no {0} is provided, this value is used directly.",
                              FieldSpecificationBase::viewKeyStruct::functionNameString() ) );
}

REGISTER_CATALOG_ENTRY( FieldSpecificationBase, SourceFluxBoundaryCondition, string const &, Group * const )

} /* namespace geos */
