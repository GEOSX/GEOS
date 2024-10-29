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
 * @file SurfaceElementRegion.hpp
 *
 */

#ifndef GEOS_MESH_SURFACEELEMENTREGION_HPP_
#define GEOS_MESH_SURFACEELEMENTREGION_HPP_

#include "ElementRegionBase.hpp"
#include "codingUtilities/EnumStrings.hpp"

namespace geos
{


class EdgeManager;

/**
 * @class SurfaceElementRegion
 *
 * The SurfaceElementRegion class contains the functionality to support the concept of a SurfaceElementRegion in the
 * element hierarchy. SurfaceElementRegion derives from ElementRegion and has an entry in the ObjectManagerBase
 * catalog.
 */
class SurfaceElementRegion : public ElementRegionBase
{
public:

  /**
   * @enum SurfaceSubRegionType
   *
   * The options for the surface subregion type
   */
  enum class SurfaceSubRegionType : integer
  {
    faceElement,
    embeddedElement
  };

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name the name of the object in the data hierarchy.
   * @param parent a pointer to the parent group in the data hierarchy.
   */
  SurfaceElementRegion( string const & name, Group * const parent );

  /**
   * @brief Deleted default constructor.
   */
  SurfaceElementRegion() = delete;

  /**
   * @brief Default destructor.
   */
  virtual ~SurfaceElementRegion() override;

  ///@}

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief Get the key name for the SurfaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static string catalogName()
  { return "SurfaceElementRegion"; }

  virtual string getCatalogName() const override final
  { return catalogName(); }

  ///@}


  /**
   * @name Generation of the face element mesh region
   */
  ///@{

  virtual void generateMesh( Group const & faceBlocks ) override;

  /**
   * @brief This function generates and adds entries to the face/fracture mesh.
   * @param time_np1 rupture time
   * @param faceManager pointer to the FaceManager object.
   * @param originalFaceToEdges face-to-edge map before the rupture.
   * @param faceIndices the local indices of the new faces that define the face element.
   * @return the local index of the new FaceElement entry.
   */
  localIndex addToFractureMesh( real64 const time_np1,
                                FaceManager const * const faceManager,
                                ArrayOfArraysView< localIndex const > const & originalFaceToEdges,
                                localIndex const faceIndices[2] );

  ///@}


  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get default aperture value.
   * @return default aperture value
   */
  real64 getDefaultAperture() const { return m_defaultAperture; }

  /**
   * @brief Get subRegion type.
   * @return subRegion type
   */
  SurfaceSubRegionType subRegionType() const { return m_subRegionType; }

  /**
   * @brief Returns the unique sub-region of type @p SUBREGION_TYPE for the current @p SurfaceElementRegion.
   * @tparam SUBREGION_TYPE The type of the sub region we're looking for.
   * @return The unique sub region.
   * @note Kills the simulation if the sub-region is not unique.
   */
  template< typename SUBREGION_TYPE >
  SUBREGION_TYPE & getUniqueSubRegion()
  {
    return getSubRegion< SUBREGION_TYPE >( getUniqueSubRegionName< SUBREGION_TYPE >() );
  }

  /**
   * @brief Returns the unique sub-region of type @p SUBREGION_TYPE for the current @p SurfaceElementRegion.
   * @tparam SUBREGION_TYPE The type of the sub region we're looking for.
   * @return The unique sub region.
   * @note Kills the simulation if the sub-region is not unique.
   */
  template< typename SUBREGION_TYPE >
  SUBREGION_TYPE const & getUniqueSubRegion() const
  {
    return getSubRegion< SUBREGION_TYPE >( getUniqueSubRegionName< SUBREGION_TYPE >() );
  }

  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    /// @return subRegion type string
    static constexpr char const * subRegionTypeString() { return "subRegionType"; }

    /// @return Fracture set string
    static constexpr char const * fractureSetString() { return "fractureSet"; }

    /// @return Default fracture aperture
    static constexpr char const * defaultApertureString() { return "defaultAperture"; }

    /// @return Rupture time string
    static constexpr char const * ruptureTimeString() { return "ruptureTime"; }

    /// @return Face block string
    static constexpr char const * faceBlockString() { return "faceBlock"; }
  };

protected:
  virtual void initializePreSubGroups() override;

private:

  /**
   * @brief Returns the name of the unique sub-region.
   * @tparam SUBREGION_TYPE The type of the sub region we're looking for.
   * @return The name unique sub region.
   * @note Kills the simulation if the sub-region is not unique.
   */
  template< typename SUBREGION_TYPE >
  string getUniqueSubRegionName() const
  {
    std::vector< string > subRegionNames;
    forElementSubRegions< SUBREGION_TYPE >( [&]( SUBREGION_TYPE const & sr )
    {
      subRegionNames.push_back( sr.getName() );
    } );
    GEOS_ERROR_IF( subRegionNames.size() != 1,
                   "Surface region \"" << getDataContext() <<
                   "\" should have one unique sub region (" << subRegionNames.size() << " found)." );
    return subRegionNames.front();
  }

  SurfaceSubRegionType m_subRegionType;

  real64 m_defaultAperture;

  /**
   * @brief One @p SurfaceElementRegion being made of one single sub-region,
   * we get the name of the corresponding @p FaceBlockABC
   * that will contain the mesh data to be imported.
   */
  string m_faceBlockName;
};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( SurfaceElementRegion::SurfaceSubRegionType,
              "faceElement",
              "embeddedElement" );

} /* namespace geos */

#endif /* CORECOMPONENTS_MESH_SurfaceElementRegion_HPP_ */
