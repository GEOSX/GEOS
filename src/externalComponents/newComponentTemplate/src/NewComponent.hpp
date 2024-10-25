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
 * @file NewComponent.hpp
 */

#ifndef COMPONENTS_NEWCOMPONENTTEMPLATE_SRC_NEWCOMPONENT_HPP_
#define COMPONENTS_NEWCOMPONENTTEMPLATE_SRC_NEWCOMPONENT_HPP_
#include "physicsSolvers/SolverBase.hpp"


namespace geos
{
namespace dataRepository
{
class Group;
}
class DomainPartition;

class NewComponent final : public SolverBase
{
public:
  NewComponent( string const & name,
                Group * const parent );
  virtual ~NewComponent() override;

  static string catalogName() { return "NewComponent"; }
  string getCatalogName() const override { return catalogName(); }

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

private:
  NewComponent() = delete;
  NewComponent( const NewComponent & ) = delete;
  NewComponent( const NewComponent && ) = delete;
  NewComponent & operator=( const NewComponent & ) = delete;
  NewComponent & operator=( const NewComponent && ) = delete;
};

} /* namespace geos */

#endif /* COMPONENTS_NEWCOMPONENTTEMPLATE_SRC_NEWCOMPONENT_HPP_ */
