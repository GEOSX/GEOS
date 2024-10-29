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
 * @file SolidMechanicsStatistics.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSTATISTICS_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSTATISTICS_HPP_

#include "physicsSolvers/FieldStatisticsBase.hpp"

namespace geos
{

class SolidMechanicsLagrangianFEM;

/**
 * @class SolidMechanicsStatistics
 *
 * Task class allowing for the computation of aggregate statistics in solid mechanics simulations
 */
class SolidMechanicsStatistics : public FieldStatisticsBase< SolidMechanicsLagrangianFEM >
{
public:

  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  SolidMechanicsStatistics( const string & name,
                            Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "SolidMechanicsStatistics"; }

  /**
   * @defgroup Tasks Interface Functions
   *
   * This function implements the interface defined by the abstract TaskBase class
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**@}*/

  /**
   * @brief Compute node-based statistics on the reservoir
   * @param[in] mesh the mesh level object
   */
  void computeNodeStatistics( MeshLevel & mesh, real64 const time ) const;

private:

  using Base = FieldStatisticsBase< SolidMechanicsLagrangianFEM >;

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    /// String for the node-based statistics
    constexpr static char const * nodeStatisticsString() { return "nodeStatistics"; }
  };

  struct NodeStatistics
  {
    /// minimum displacement in each direction (x, y, z)
    array1d< real64 > minDisplacement;
    /// maximum displacement in each direction (x, y, z)
    array1d< real64 > maxDisplacement;
  };

  void registerDataOnMesh( Group & meshBodies ) override;

};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSTATISTICS_HPP_ */
