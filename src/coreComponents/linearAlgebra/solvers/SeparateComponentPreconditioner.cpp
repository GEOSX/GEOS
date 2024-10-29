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
 * @file SeparateComponentPreconditioner.cpp
 */

#include "SeparateComponentPreconditioner.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

namespace geos
{

template< typename LAI >
SeparateComponentPreconditioner< LAI >::
SeparateComponentPreconditioner( localIndex const numComp,
                                 std::unique_ptr< PreconditionerBase< LAI > > precond )
  : Base(),
  m_numComp( numComp ),
  m_precond( std::move( precond ) )
{
  GEOS_LAI_ASSERT_GT( m_numComp, 0 );
  GEOS_LAI_ASSERT( m_precond );
}

template< typename LAI >
SeparateComponentPreconditioner< LAI >::~SeparateComponentPreconditioner() = default;

template< typename LAI >
void SeparateComponentPreconditioner< LAI >::setup( Matrix const & mat )
{
  Base::setup( mat );
  // TODO: if matrix structure hasn't changed, can just copy entries into existing m_matSC
  mat.separateComponentFilter( m_matSC, m_numComp );
  m_precond->setup( m_matSC );
}

template< typename LAI >
void SeparateComponentPreconditioner< LAI >::apply( Vector const & src,
                                                    Vector & dst ) const
{
  m_precond->apply( src, dst );
}

template< typename LAI >
void SeparateComponentPreconditioner< LAI >::clear()
{
  Base::clear();
  m_precond->clear();
  m_matSC.reset();
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOS_USE_TRILINOS
template class SeparateComponentPreconditioner< TrilinosInterface >;
#endif

#ifdef GEOS_USE_HYPRE
template class SeparateComponentPreconditioner< HypreInterface >;
#endif

#ifdef GEOS_USE_PETSC
template class SeparateComponentPreconditioner< PetscInterface >;
#endif

}
