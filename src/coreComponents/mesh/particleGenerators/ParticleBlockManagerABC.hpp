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

#ifndef GEOS_PARTICLEBLOCKMANAGERABC_HPP
#define GEOS_PARTICLEBLOCKMANAGERABC_HPP

#include "dataRepository/Group.hpp"
#include "mesh/ParticleType.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"

namespace geos
{

/**
 * @brief Abstract base class for ParticleBlockManager.
 */
class ParticleBlockManagerABC : public dataRepository::Group
{
public:

  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param parent The parent Group.
   */
  ParticleBlockManagerABC( string const & name, Group * const parent ):
    Group( name, parent )
  {
    // Left blank
  }

  /**
   * @brief Returns a group containing the cell blocks as ParticleBlockABC instances
   * @return Mutable reference to the cell blocks group.
   *
   * @note It should probably be better not to expose a non-const accessor here.
   */
  virtual Group & getParticleBlocks() = 0;

  /**
   * @brief Returns a group containing the cell blocks as ParticleBlockABC instances
   * @return Const reference to the Group instance.
   */
  virtual const Group & getParticleBlocks() const = 0;

};

}
#endif // include guard
