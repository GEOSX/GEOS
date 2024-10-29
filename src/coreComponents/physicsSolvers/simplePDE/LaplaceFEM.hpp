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

#ifndef GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_FEM_HPP_
#define GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_FEM_HPP_

#include "physicsSolvers/simplePDE/LaplaceBaseH1.hpp"  // a base class shared by all Laplace solvers

namespace geos
{

// Like most physics solvers, the Laplace solver derives from a generic SolverBase class.
// The base class is densely Doxygen-commented and worth a look if you have not done so already.
// Most important system assembly steps, linear and non-linear resolutions, and time-stepping mechanisms
// are implemented at the SolverBase class level and can thus be used in Laplace without needing reimplementation.

//START_SPHINX_INCLUDE_BEGINCLASS
class LaplaceFEM : public LaplaceBaseH1
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  LaplaceFEM() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the
  /// tree structure of classes)
  LaplaceFEM( const string & name,
              Group * const parent );

  /// Destructor
  virtual ~LaplaceFEM() override;

  /// "CatalogName()" return the string used as XML tag in the input file.  It ties the XML tag with
  /// this C++ classes. This is important.
  static string catalogName() { return "LaplaceFEM"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

//END_SPHINX_INCLUDE_BEGINCLASS
// /**
//  * @defgroup Solver Interface Functions
//  *
//  * These functions provide the primary interface that is required for derived classes
//  */
// /**@{*/

//START_SPHINX_INCLUDE_SOLVERINTERFACE
  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = false ) override;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

//END_SPHINX_INCLUDE_SOLVERINTERFACE

  /**@}*/

};
} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_FEM_HPP_ */
