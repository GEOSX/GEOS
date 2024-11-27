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
 * @file FaultTractionUpdateBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_FAULTTRACTIONUPDATE_BASE_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_FAULTTRACTIONUPDATE_BASE_HPP

#include "common/DataTypes.hpp"

// #include "mesh/SurfaceElementRegion.hpp"

namespace geos
{

class DomainPartition;
class SurfaceElementSubRegion;
namespace inducedSeismicity
{
class FaultTractionUpdateBase
{
public:

  virtual ~FaultTractionUpdateBase() = default;

  virtual real64 updateFaultTraction( real64 const & time_n,
                                      real64 const & dt,
                                      const int cycleNumber,
                                      DomainPartition & domain ) const

  {
    GEOS_UNUSED_VAR( time_n, cycleNumber, domain );
    return dt;
  }

  virtual void registerMissingDataOnMesh( SurfaceElementSubRegion & subRegion, string const & solverName ) const
  {
    GEOS_UNUSED_VAR( subRegion, solverName );
  }

};

enum class TractionUpdateType : integer
{
  SpringSlider,
  OneWayCoupled
};


} // namespace inducedSeismicity

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_STRESS_SOLVER_WRAPPER_BASE_HPP */
