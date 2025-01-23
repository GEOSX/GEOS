/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#ifndef GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEARTHQUAKE_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEARTHQUAKE_HPP

#include "ImplicitQDRateAndState.hpp"

namespace geos
{

template< typename RSSOLVER_TYPE = ImplicitQDRateAndState >
class QuasiDynamicEarthQuake : public RSSOLVER_TYPE
{
public:

  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  QuasiDynamicEarthQuake() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  QuasiDynamicEarthQuake( const string & name,
                          dataRepository::Group * const parent );

  /// Destructor
  virtual ~QuasiDynamicEarthQuake() override;

  static string catalogName() { return RSSOLVER_TYPE::derivedSolverPrefix() + "QuasiDynamicEQ"; }

  /**
   * @return Get the final class Catalog name
   */
  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public RSSOLVER_TYPE::viewKeyStruct
  {
    /// stress solver name
    static constexpr char const * stressSolverNameString() { return "stressSolverName"; }
  };

  void postInputInitialization() override final;

  void setTargetDispJump( DomainPartition & domain ) const;

  virtual real64 updateStresses( real64 const & time_n,
                                 real64 const & dt,
                                 const int cycleNumber,
                                 DomainPartition & domain ) const override final;

private:

  string m_stressSolverName;

  PhysicsSolverBase * m_stressSolver;

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_QUASIDYNAMICEQBASE_HPP */
