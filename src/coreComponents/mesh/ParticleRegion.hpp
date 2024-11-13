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
 * @file ParticleRegion.hpp
 *
 */

#ifndef GEOS_MESH_PARTICLEREGION_HPP_
#define GEOS_MESH_PARTICLEREGION_HPP_

#include "ParticleRegionBase.hpp"

namespace geos
{
/**
 * @class ParticleRegion
 *
 * The ParticleRegion class contains the functionality to support the concept of a ParticleRegion in the element
 * hierarchy. ParticleRegion derives from ParticleRegionBase and has an entry in the ObjectManagerBase catalog.
 *
 *
 */
class ParticleRegion : public ParticleRegionBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{


  /**
   * @brief Constructor.
   * @param name the name of the object in the data hierarchy.
   * @param parent a pointer to the parent group in the data hierarchy.
   */
  ParticleRegion( string const & name, Group * const parent );

  /**
   * @brief Deleted default constructor.
   */
  ParticleRegion() = delete;

  /**
   * @brief Destructor.
   */
  virtual ~ParticleRegion() override;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief The key name for the FaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static string catalogName()
  { return "ParticleRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual string getCatalogName() const override final
  { return ParticleRegion::catalogName(); }

  ///@}

  /**
   * @name Generation of the particle region
   */
  ///@{

  /**
   * @brief Add a particleBlockRegion name to the list.
   * @param particleBlockName string containing the cell block region name.
   */
  void addParticleBlockName( string const & particleBlockName )
  {
    m_particleBlockNames.emplace_back( particleBlockName );
  }

  /**
   * @brief Get the list of particleBlock names.
   * @return A string array of particle block names
   */
  string_array getParticleBlockNames()
  {
    return m_particleBlockNames;
  }

  virtual void generateMesh( Group & particleBlocks ) override;

  ///@}

  /**
   * @brief Calculate and return the locations of all particle corners
   * @return The list of particle corner coordinates
   */
  array2d< real64 > getParticleCorners() const;

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ParticleRegionBase::viewKeyStruct
  {
    /// @return String key for the cell block names
    static constexpr char const * sourceParticleBlockNamesString() {return "particleBlocks"; }
  };


private:

  // Cell block names
  string_array m_particleBlockNames;

};

} /* namespace geos */

#endif /* GEOS_MESH_PARTICLEREGION_HPP_ */
