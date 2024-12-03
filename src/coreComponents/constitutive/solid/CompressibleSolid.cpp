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
 * @file CompressibleSolid.cpp
 */

#include "CompressibleSolid.hpp"
#include "porosity/PressurePorosity.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"
#include "constitutive/permeability/CarmanKozenyPermeability.hpp"
#include "constitutive/permeability/ExponentialDecayPermeability.hpp"
#include "constitutive/permeability/ParallelPlatesPermeability.hpp"
#include "constitutive/permeability/PressurePermeability.hpp"
#include "constitutive/permeability/SlipDependentPermeability.hpp"
#include "constitutive/permeability/WillisRichardsPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

template< typename PORO_TYPE,
          typename PERM_TYPE >
CompressibleSolid< PORO_TYPE, PERM_TYPE >::CompressibleSolid( string const & name, Group * const parent ):
  CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >( name, parent )
{}

template< typename PORO_TYPE,
          typename PERM_TYPE >
CompressibleSolid< PORO_TYPE, PERM_TYPE >::~CompressibleSolid() = default;

// Register all CompressibleSolid model types.
typedef CompressibleSolid< PressurePorosity, ConstantPermeability > CompressibleRockConstant;
typedef CompressibleSolid< PressurePorosity, CarmanKozenyPermeability > CompressibleRockCK;
typedef CompressibleSolid< PressurePorosity, PressurePermeability > CompressibleRockPressurePerm;
typedef CompressibleSolid< PressurePorosity, ExponentialDecayPermeability > FaultED;
typedef CompressibleSolid< PressurePorosity, ParallelPlatesPermeability > FractureRock;
typedef CompressibleSolid< PressurePorosity, SlipDependentPermeability > Fault;
typedef CompressibleSolid< PressurePorosity, WillisRichardsPermeability > FaultWR;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRockConstant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRockCK, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRockPressurePerm, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FractureRock, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FaultED, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, Fault, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, FaultWR, string const &, Group * const )


}
} /* namespace geos */
