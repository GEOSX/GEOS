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
 * @file ElementRegionManager.hpp
 */

#ifndef GEOS_MESH_ELEMENTREGIONMANAGER_HPP
#define GEOS_MESH_ELEMENTREGIONMANAGER_HPP

#include "constitutive/ConstitutiveManager.hpp"
#include "CellElementRegion.hpp"
#include "CellElementSubRegion.hpp"
#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "dataRepository/ReferenceWrapper.hpp"
#include "SurfaceElementRegion.hpp"
#include "WellElementRegion.hpp"

#include "common/TypeDispatch.hpp"

namespace geos
{

class MeshManager;

/**
 * @class ElementRegionManager
 * @brief The ElementRegionManager class provides an interface to ObjectManagerBase in order to manage ElementRegion
 * data
 */
class ElementRegionManager : public ObjectManagerBase
{
public:

  /**
   * @brief Group key associated with elementRegionsGroup.
   */
  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    /// @return element regions group string key
    static constexpr auto elementRegionsGroup() { return "elementRegionsGroup"; }
  };

  /**
   * Limit on max number of nodes for each element
   */
  constexpr static int maxNumNodesPerElem = 8;

  /**
   * @brief The ElementViewAccessor at the ElementRegionManager level is an array of array of VIEWTYPE.
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ElementViewAccessor = array1d< array1d< VIEWTYPE > >;

  /**
   * @brief The ElementViewAccessor at the ElementRegionManager level is the
   *   type resulting from ElementViewAccessor< VIEWTYPE >::toNestedView().
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ElementView = typename ElementViewAccessor< VIEWTYPE >::NestedViewType;

  /**
   * @brief The ElementViewAccessor at the ElementRegionManager level is the
   *   type resulting from ElementViewAccessor< VIEWTYPE >::toNestedViewConst().
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ElementViewConst = typename ElementViewAccessor< VIEWTYPE >::NestedViewTypeConst;

  /**
   * @brief The ElementViewAccessor at the ElementRegionManager level is a 2D array of ReferenceWrapper around VIEWTYPE.
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ElementReferenceAccessor = array1d< array1d< ReferenceWrapper< VIEWTYPE > > >;

  /**
   * @brief The MaterialViewAccessor at the ElementRegionManager level is a 3D array of VIEWTYPE.
   * @tparam VIEWTYPE data type
   * var[elementRegionIndex][elementSubRegionIndex][materialIndexInRegion]
   */
  template< typename VIEWTYPE >
  using MaterialViewAccessor = array1d< array1d< array1d< VIEWTYPE > > >;

  /**
   * @brief The ConstitutiveRelationAccessor at the ElementRegionManager level is a 3D array of CONSTITUTIVE_TYPE
   * @tparam CONSTITUTIVE_TYPE constitutive type
   */
  template< typename CONSTITUTIVE_TYPE >
  using ConstitutiveRelationAccessor = array1d< array1d< array1d< CONSTITUTIVE_TYPE * > > >;

  /**
   * @brief The function is to return the name of the ElementRegionManager in the object catalog
   * @return string that contains the catalog name used to register/lookup this class in  the object catalog
   */
  static string catalogName()
  { return "ZoneManager"; }

  /**
   * @brief Virtual access to catalogName()
   * @return string that contains the catalog name used to register/lookup this class in the object catalog
   */
  virtual string getCatalogName() const override final
  { return catalogName(); }

  /**
   * @brief Constructor.
   * @param [in] name the name of this ObjectManager
   * @param [in] parent the parent Group
   */
  ElementRegionManager( string const & name, Group * const parent );

  /**
   * @brief Destructor
   */
  virtual ~ElementRegionManager() override;

  /**
   * @brief Get the number of elements in all ElementSubRegions of type T.
   * @return number of elements
   */
  template< typename T = ElementSubRegionBase >
  localIndex getNumberOfElements() const
  {
    localIndex numElem = 0;
    this->forElementSubRegions< T >( [&]( ElementSubRegionBase const & elementSubRegion )
    {
      numElem += elementSubRegion.size();
    } );
    return numElem;
  }

  /**
   * @brief Generate the mesh. Produce an error if a required cellBlock doesn't exist in the source mesh.
   * @param [in,out] cellBlockManager Reference to the abstract cell block manager.
   */
  void generateMesh( CellBlockManagerABC const & cellBlockManager );

  /**
   * @brief Generate the wells.
   * @param [in] cellBlockManager pointer to cellBlockManager
   * @param [in] meshLevel pointer to meshLevel
   */
  void generateWells( CellBlockManagerABC const & cellBlockManager, MeshLevel & meshLevel );

  /**
   * @brief Build sets from the node sets
   * @param[in] nodeManager The node manager that will provide the node sets.
   * @note ElementRegionManager's sub-regions need to be properly defined.
   */
  void buildSets( NodeManager const & nodeManager );

  /**
   * @brief Create a new ElementRegion object as a child of this group.
   * @param childKey catalog key of the new ElementRegion derived type to create
   * @param childName name of the new ElementRegion object
   * @return pointer to the created ElementRegion object
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;
//  virtual void ReadXMLsub( xmlWrapper::xmlNode const & targetNode ) override;

  /**
   * @brief Expand any catalogs in the data structure
   */
  virtual void expandObjectCatalogs() override;

  /**
   * @brief Inform the schema generator of any deviations between the xml and GEOS data structures.
   * @param schemaRoot        XML node corresponding to the root
   * @param schemaParent      XML node for the parent node
   * @param documentationType type of XML schema generated
   */
  virtual void setSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                    xmlWrapper::xmlNode schemaParent,
                                    integer documentationType ) override;

  using Group::resize;

  /**
   * @brief Set the number of elements for a set of element regions.
   * @param numElements list of the new element numbers
   * @param regionNames list of the element region names
   * @param elementTypes list of the element types
   */
  void resize( integer_array const & numElements,
               string_array const & regionNames,
               string_array const & elementTypes );

  /**
   * @brief Set the maximum local and global index.
   */
  virtual void setMaxGlobalIndex() override final;

  /**
   * @brief Get a collection of element regions
   * @return reference to immutable subGroupMap
   */
  subGroupMap const & getRegions() const
  {
    return this->getGroup( groupKeyStruct::elementRegionsGroup() ).getSubGroups();
  }

  /**
   * @brief Get a collection of element regions.
   * @return reference to mutable subGroupMap
   */
  subGroupMap & getRegions()
  {
    return this->getGroup( groupKeyStruct::elementRegionsGroup() ).getSubGroups();
  }

  /**
   * @brief Get a element region.
   * @param key The key of element region, either name or number.
   * @return Reference to const T.
   * @throw std::domain_error if the requested region doesn't exist.
   */
  template< typename T=ElementRegionBase, typename KEY_TYPE=void >
  T const & getRegion( KEY_TYPE const & key ) const
  {
    return this->getGroup( groupKeyStruct::elementRegionsGroup() ).getGroup< T >( key );
  }

  /**
   * @brief Get a element region.
   * @param key The key of the element region, either name or number.
   * @return Reference to T.
   * @throw std::domain_error if the requested region doesn't exist.
   */
  template< typename T=ElementRegionBase, typename KEY_TYPE=void >
  T & getRegion( KEY_TYPE const & key )
  {
    return this->getGroup( groupKeyStruct::elementRegionsGroup() ).getGroup< T >( key );
  }

  /**
   * @brief Determines if an ElementRegion with the input name exists.
   * @tparam T The type of ElementRegion. May be a specific derived type of ElementRegionBase.
   * @param name The name/key of the ElementRegion
   * @return true if the region exists, false if not.
   */
  template< typename T=ElementRegionBase >
  bool hasRegion( string const & name ) const
  {
    return this->getGroup( groupKeyStruct::elementRegionsGroup() ).hasGroup< T >( name );
  }

  /**
   * @brief Get number of the regions.
   * @return number of the regions
   */
  localIndex numRegions() const
  {
    return this->getRegions().size();
  }

  /**
   * @brief Produce a map from cell block indices to element region and subregion indices
   * @param cellBlockManager the CellBlocKManager
   * @return a (numBlock x 2) array with each row corresponding to a cell block and containing
   *         region (first entry) and subregion (second entry) indices, or -1 if block was not used.
   */
  array2d< localIndex > getCellBlockToSubRegionMap( CellBlockManagerABC const & cellBlockManager ) const;

  /**
   * @brief This function is used to launch kernel function over all the element regions with region type =
   * ElementRegionBase.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ElementRegionBase, typename ... REGIONTYPES, typename LAMBDA >
  void forElementRegions( LAMBDA && lambda )
  {
    this->getGroup( groupKeyStruct::elementRegionsGroup() ).forSubGroups< REGIONTYPE, REGIONTYPES... >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over all the element regions with region type =
   * ElementRegionBase.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ElementRegionBase, typename ... REGIONTYPES, typename LAMBDA >
  void forElementRegions( LAMBDA && lambda ) const
  {
    this->getGroup( groupKeyStruct::elementRegionsGroup() ).forSubGroups< REGIONTYPE, REGIONTYPES... >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the target element regions with region type =
   * ElementRegionBase.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ElementRegionBase, typename ... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    this->getGroup( groupKeyStruct::elementRegionsGroup() ).forSubGroups< REGIONTYPE, REGIONTYPES... >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the target element regions with region type =
   * ElementRegionBase.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ElementRegionBase, typename ... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    this->getGroup( groupKeyStruct::elementRegionsGroup() ).forSubGroups< REGIONTYPE, REGIONTYPES... >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over all the types of element regions.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forElementRegionsComplete( LAMBDA lambda ) const
  {
    forElementRegionsComplete< CellElementRegion, SurfaceElementRegion,
                               WellElementRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over all the types of element regions.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forElementRegionsComplete( LAMBDA lambda )
  {
    forElementRegionsComplete< CellElementRegion, SurfaceElementRegion,
                               WellElementRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over all the element regions that can be casted to one of
   * the specified region types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename REGIONTYPE, typename ... REGIONTYPES, typename LAMBDA >
  void forElementRegionsComplete( LAMBDA lambda )
  {
    for( localIndex er=0; er<this->numRegions(); ++er )
    {
      ElementRegionBase & elementRegion = this->getRegion( er );

      Group::applyLambdaToContainer< REGIONTYPE, REGIONTYPES... >( elementRegion, [&]( auto & castedRegion )
      {
        lambda( er, castedRegion );
      } );
    }
  }

  /**
   * @brief This const function is used to launch kernel function over all the element regions that can be casted to one
   * of the specified region types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename REGIONTYPE, typename ... REGIONTYPES, typename LAMBDA >
  void forElementRegionsComplete( LAMBDA lambda ) const
  {
    for( localIndex er=0; er<this->numRegions(); ++er )
    {
      ElementRegionBase const & elementRegion = this->getRegion( er );

      Group::applyLambdaToContainer< REGIONTYPE, REGIONTYPES... >( elementRegion, [&]( auto const & castedRegion )
      {
        lambda( er, castedRegion );
      } );
    }
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element regions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda ) const
  {
    forElementRegionsComplete< CellElementRegion, SurfaceElementRegion,
                               WellElementRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element regions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda )
  {
    forElementRegionsComplete< CellElementRegion, SurfaceElementRegion,
                               WellElementRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element regions with region type =
   * specified element region types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename REGIONTYPE, typename ... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda )
  {
    forElementRegions< REGIONTYPE, REGIONTYPES... >( targetRegions, [&] ( localIndex const targetIndex,
                                                                          auto & elementRegion )
    {
      lambda( targetIndex, elementRegion.getIndexInParent(), elementRegion );
    } );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element regions with region
   * type = specified element region types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename REGIONTYPE, typename ... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda ) const
  {
    forElementRegions< REGIONTYPE, REGIONTYPES... >( targetRegions, [&] ( localIndex const targetIndex,
                                                                          auto const & elementRegion )
    {
      lambda( targetIndex, elementRegion.getIndexInParent(), elementRegion );
    } );
  }

  /**
   * @brief This function is used to launch kernel function over the element subregions of all the subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    forElementSubRegions< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion,
                          WellElementSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the element subregions of all the subregion
   * types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda ) const
  {
    forElementSubRegions< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion,
                          WellElementSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element subregions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forElementSubRegions< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion,
                          WellElementSubRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element subregions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forElementSubRegions< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion,
                          WellElementSubRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the element subregions of the specified subregion
   * types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda )
  {
    forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >(
      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const,
                                                   localIndex const,
                                                   ElementRegionBase &,
                                                   auto & subRegion )
    {
      lambda( subRegion );
    }
      );
  }

  /**
   * @brief This const function is used to launch kernel function over the element subregions of the specified subregion
   * types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegions( LAMBDA && lambda ) const
  {
    forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >(
      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const,
                                                   localIndex const,
                                                   ElementRegionBase const &,
                                                   auto const & subRegion )
    {
      lambda( subRegion );
    } );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element subregions with the
   * specified subregion types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >( targetRegions,
                                                                      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const targetIndex,
                                                                                                                   localIndex const,
                                                                                                                   localIndex const,
                                                                                                                   ElementRegionBase &,
                                                                                                                   auto & subRegion )
    {
      lambda( targetIndex, subRegion );
    } );
  }

  template< typename ... SUBREGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegions( types::TypeList< SUBREGIONTYPES... >, LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forElementSubRegionsComplete< SUBREGIONTYPES... >( targetRegions,
                                                                      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const targetIndex,
                                                                                                                   localIndex const,
                                                                                                                   localIndex const,
                                                                                                                   ElementRegionBase &,
                                                                                                                   auto & subRegion )
    {
      lambda( targetIndex, subRegion );
    } );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element subregions with the
   * specified subregion types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forElementSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >( targetRegions,
                                                                      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const targetIndex,
                                                                                                                   localIndex const,
                                                                                                                   localIndex const,
                                                                                                                   ElementRegionBase const &,
                                                                                                                   auto const & subRegion )
    {
      lambda( targetIndex, subRegion );
    } );
  }

  /**
   * @brief This const function is used to launch kernel function over the element subregions of all subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forElementSubRegionsComplete( LAMBDA && lambda ) const
  {
    forElementSubRegionsComplete< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion,
                                  WellElementSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the element subregions of all subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forElementSubRegionsComplete( LAMBDA && lambda )
  {
    forElementSubRegionsComplete< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion,
                                  WellElementSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element subregions
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forElementSubRegionsComplete< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion, WellElementSubRegion >( targetRegions,
                                                                                                                                std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element subregions
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forElementSubRegionsComplete< CellElementSubRegion, FaceElementSubRegion, EmbeddedSurfaceSubRegion, WellElementSubRegion >( targetRegions,
                                                                                                                                std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over all the element subregions that can be casted to one of
   * the specified subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegionsComplete( LAMBDA && lambda )
  {
    for( localIndex er=0; er<this->numRegions(); ++er )
    {
      ElementRegionBase & elementRegion = this->getRegion( er );

      for( localIndex esr=0; esr<elementRegion.numSubRegions(); ++esr )
      {
        ElementSubRegionBase & subRegion = elementRegion.getSubRegion( esr );

        Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto & castedSubRegion )
        {
          lambda( er, esr, elementRegion, castedSubRegion );
        } );
      }
    }
  }

  /**
   * @brief This const function is used to launch kernel function over all the element subregions that can be casted to
   * one of the specified subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forElementSubRegionsComplete( LAMBDA && lambda ) const
  {
    for( localIndex er=0; er<this->numRegions(); ++er )
    {
      ElementRegionBase const & elementRegion = this->getRegion( er );

      for( localIndex esr=0; esr<elementRegion.numSubRegions(); ++esr )
      {
        ElementSubRegionBase const & subRegion = elementRegion.getSubRegion( esr );

        Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto const & castedSubRegion )
        {
          lambda( er, esr, elementRegion, castedSubRegion );
        } );
      }
    }
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element subregions that can be
   * casted to one of the specified subregion types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forElementRegions( targetRegions, [&] ( localIndex const targetIndex, ElementRegionBase & elementRegion )
    {
      localIndex const er = elementRegion.getIndexInParent();

      if( er>-1 )
      {
        for( localIndex esr=0; esr<elementRegion.numSubRegions(); ++esr )
        {
          ElementSubRegionBase & subRegion = elementRegion.getSubRegion( esr );

          Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto & castedSubRegion )
          {
            lambda( targetIndex, er, esr, elementRegion, castedSubRegion );
          } );
        }
      }
    } );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element subregions that can
   * be casted to one of the specified subregion types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forElementSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forElementRegions( targetRegions, [&] ( localIndex const targetIndex, ElementRegionBase const & elementRegion )
    {
      localIndex const er = elementRegion.getIndexInParent();

      if( er>-1 )
      {
        for( localIndex esr=0; esr<elementRegion.numSubRegions(); ++esr )
        {
          ElementSubRegionBase const & subRegion = elementRegion.getSubRegion( esr );

          Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto const & castedSubRegion )
          {
            lambda( targetIndex, er, esr, elementRegion, castedSubRegion );
          } );
        }
      }
    } );
  }


  /**
   * @brief This is a const function to construct a ElementViewAccessor to access the data registered on the mesh.
   * @tparam FIELD_TRAIT field type
   * @param neighborName neighbor data name
   * @return ElementViewAccessor that contains traits::ViewTypeConst< typename TRAIT::type > data
   */
  template< typename FIELD_TRAIT >
  ElementViewAccessor< traits::ViewTypeConst< typename FIELD_TRAIT::type > >
  constructFieldAccessor( string const & neighborName = string() ) const;

  /**
   * @brief This is a const function to construct a ElementViewAccessor to access the data registered on the mesh.
   * @tparam VIEWTYPE data type
   * @param name view name of the data
   * @param neighborName neighbor data name
   * @return ElementViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ElementViewAccessor< LHS >
  constructViewAccessor( string const & name, string const & neighborName = string() ) const;

  /**
   * @brief This is a function to construct a ElementViewAccessor to access the data registered on the mesh.
   * @tparam VIEWTYPE data type
   * @param name view name of the data
   * @param neighborName neighbor data name
   * @return ElementViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ElementViewAccessor< LHS >
  constructViewAccessor( string const & name, string const & neighborName = string() );

  /**
   * @brief This is a function to construct a ElementViewAccessor to access array data registered on the mesh.
   * @tparam T data type
   * @tparam NDIM number of array dimensions
   * @tparam PERM layout permutation sequence type
   * @param name view name of the data
   * @param neighborName neighbor data name
   * @return ElementViewAccessor that contains ArrayView<T const, NDIM> of data
   */
  template< typename T, int NDIM, typename PERM = defaultLayout< NDIM > >
  ElementViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
  constructArrayViewAccessor( string const & name, string const & neighborName = string() ) const;

  /**
   * @brief This is a const function to construct a ElementViewAccessor to access the data registered on the mesh.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param neighborName neighbor data name
   * @return ElementViewAccessor that contains pointers to wrapped VIEWTYPE data
   */
  template< typename VIEWTYPE >
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
  constructReferenceAccessor( string const & viewName, string const & neighborName = string() ) const;

  /**
   * @brief This is a function to construct a ElementViewAccessor to access the data registered on the mesh.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param neighborName neighbor data name
   * @return ElementViewAccessor that contains pointers to wrapped VIEWTYPE data
   */
  template< typename VIEWTYPE >
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
  constructReferenceAccessor( string const & viewName, string const & neighborName = string() );

  /**
   * @brief This is a const function to construct a MaterialViewAccessor to access the material data.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param cm pointer to ConstitutiveManager
   * @return MaterialViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  MaterialViewAccessor< LHS >
  constructFullMaterialViewAccessor( string const & viewName,
                                     constitutive::ConstitutiveManager const & cm ) const;

  /**
   * @brief This is a function to construct a MaterialViewAccessor to access the material data.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param cm pointer to ConstitutiveManager
   * @return MaterialViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  MaterialViewAccessor< LHS >
  constructFullMaterialViewAccessor( string const & viewName,
                                     constitutive::ConstitutiveManager const & cm );

  /**
   * @brief This is a const function to construct a MaterialViewAccessor to access the material data for specified
   * regions/materials.
   * @tparam FIELD_TRAIT field trait
   * @param regionNames list of region names
   * @param materialNames list of corresponding material names
   * @param allowMissingViews flag to indicate whether it is allowed to miss the specified material data in material
   * list
   * @return ElementViewAccessor that contains traits::ViewTypeConst< typename FIELD_TRAIT::type > data
   */
  template< typename FIELD_TRAIT >
  ElementViewAccessor< traits::ViewTypeConst< typename FIELD_TRAIT::type > >
  constructMaterialFieldAccessor( arrayView1d< string const > const & regionNames,
                                  arrayView1d< string const > const & materialNames,
                                  bool const allowMissingViews = false ) const;

  /**
   * @brief This is a const function to construct a MaterialViewAccessor to access the material data for specified
   * material type.
   * @tparam MATERIAL_TYPE base type of material model
   * @tparam FIELD_TRAIT field trait
   * @param allowMissingViews flag to indicate whether it is allowed to miss the specified material data in material
   * list
   * @return ElementViewAccessor that contains traits::ViewTypeConst< typename TRAIT::type > data
   */
  template< typename MATERIAL_TYPE, typename FIELD_TRAIT >
  ElementViewAccessor< traits::ViewTypeConst< typename FIELD_TRAIT::type > >
  constructMaterialFieldAccessor( bool const allowMissingViews = false ) const;


  /**
   * @brief This is a const function to construct a MaterialViewAccessor to access the material data for specified
   * regions/materials.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param regionNames list of region names
   * @param materialKeyName key of the wrapper that holds the material name on the subRegion.
   * @param allowMissingViews flag to indicate whether it is allowed to miss the specified material data in material
   * list
   * @return ElementViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ElementViewAccessor< LHS >
  constructMaterialViewAccessor( string const & viewName,
                                 arrayView1d< string const > const & regionNames,
                                 string const & materialKeyName,
                                 bool const allowMissingViews = false ) const;

  /**
   * @brief This is a function to construct a MaterialViewAccessor to access the material data for specified
   * regions/materials.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param regionNames list of region names
   * @param materialKeyName key of the wrapper that holds the material name on the subRegion.
   * @param allowMissingViews flag to indicate whether it is allowed to miss the specified material data in material
   * list
   * @return ElementViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ElementViewAccessor< LHS >
  constructMaterialViewAccessor( string const & viewName,
                                 arrayView1d< string const > const & regionNames,
                                 string const & materialKeyName,
                                 bool const allowMissingViews = false );

  /**
   * @brief Construct a view accessor for material data, assuming array as storage type
   * @tparam T underlying data type
   * @tparam NDIM number of array dimensions
   * @tparam PERM layout permutation sequence type
   * @param viewName view name of the data
   * @param regionNames list of region names
   * @param materialKeyName key of the wrapper that holds the material name on the subRegion.
   * @param allowMissingViews flag to indicate whether it is allowed to miss the specified material data in material list
   * @return MaterialViewAccessor that contains the data views
   */
  template< typename T, int NDIM, typename PERM = defaultLayout< NDIM > >
  ElementViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
  constructMaterialArrayViewAccessor( string const & viewName,
                                      arrayView1d< string const > const & regionNames,
                                      string const & materialKeyName,
                                      bool const allowMissingViews = false ) const;

  /**
   * @brief Construct a const view accessor to material data for specified material type.
   * @tparam MATERIALTYPE base type of material model
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @return ElementViewAccessor that contains VIEWTYPE data. Empty views are returned
   *         for subregions that don't contain a model derived from MODELTYPE.
   */
  template< typename MATERIALTYPE, typename VIEWTYPE, typename LHS=VIEWTYPE >
  ElementViewAccessor< LHS >
  constructMaterialViewAccessor( string const & viewName ) const;

  /**
   * @brief Construct a const view accessor for material data, assuming array as storage type
   * @tparam MATERIALTYPE
   * @tparam T underlying data type
   * @tparam NDIM number of array dimensions
   * @tparam PERM layout permutation sequence type
   * @param viewName view name of the data
   * @return MaterialViewAccessor that contains the data views. Empty views are returned
   *         for subregions that don't contain a model derived from MODELTYPE.
   */
  template< typename MATERIALTYPE, typename T, int NDIM, typename PERM = defaultLayout< NDIM > >
  ElementViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
  constructMaterialArrayViewAccessor( string const & viewName ) const;

  /**
   * @brief Construct a ConstitutiveRelationAccessor.
   * @tparam CONSTITUTIVE_TYPE constitutive type
   * @param cm pointer to ConstitutiveManager
   * @return ConstitutiveRelationAccessor
   */
  template< typename CONSTITUTIVE_TYPE >
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
  constructFullConstitutiveAccessor( constitutive::ConstitutiveManager const & cm ) const;


  /**
   * @brief Construct a ConstitutiveRelationAccessor.
   * @tparam CONSTITUTIVE_TYPE constitutive type
   * @param cm pointer to ConstitutiveManager
   * @return ConstitutiveRelationAccessor
   */
  template< typename CONSTITUTIVE_TYPE >
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
  constructFullConstitutiveAccessor( constitutive::ConstitutiveManager const & cm );

  using Group::packSize;
  using Group::pack;
  using ObjectManagerBase::packGlobalMapsSize;
  using ObjectManagerBase::packGlobalMaps;
  using ObjectManagerBase::unpackGlobalMaps;
  using ObjectManagerBase::packUpDownMapsSize;
  using ObjectManagerBase::packUpDownMaps;
  using ObjectManagerBase::unpackUpDownMaps;

  /**
   * @brief Get the buffer size needed to pack all the wrappers of all the sub regions of all the regions.
   * @param packList list of indices to pack
   * @return the size of the buffer required to pack the wrappers
   */
  int packSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack all the wrappers of all the sub regions of all the regions.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of data packed to the buffer
   */
  int pack( buffer_unit_type * & buffer,
            ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /// @copydoc dataRepository::Group::unpack
  using ObjectManagerBase::unpack;

  /**
   * @brief Unpack a buffer.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to unpack
   * @return the size of data unpacked
   */
  int unpack( buffer_unit_type const * & buffer,
              ElementViewAccessor< arrayView1d< localIndex > > & packList );

  /**
   * @brief Unpack a buffer.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to unpack
   * @return the size of data unpacked.
   */
  int unpack( buffer_unit_type const * & buffer,
              ElementReferenceAccessor< array1d< localIndex > > & packList );

  /**
   * @brief Get the size of the buffer to be packed.
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  int packGlobalMapsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack a buffer.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  int packGlobalMaps( buffer_unit_type * & buffer,
                      ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Unpack a buffer.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to pack
   * @return the size of the data unpacked
   */
  int unpackGlobalMaps( buffer_unit_type const * & buffer,
                        ElementViewAccessor< ReferenceWrapper< localIndex_array > > & packList );

  /**
   * @brief Get the buffer size needed to pack element-to-node and element-to-face maps.
   * @param packList list of indices to pack
   * @return the size of data packed.
   */
  int packUpDownMapsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Get the buffer size needed to pack element-to-node and element-to-face maps.
   * @param packList list of indices to pack
   * @return the size of data packed.
   */
  int packUpDownMapsSize( ElementReferenceAccessor< array1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack element-to-node and element-to-face maps.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of data packed.
   */
  int packUpDownMaps( buffer_unit_type * & buffer,
                      ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack element-to-node and element-to-face maps.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of data packed.
   */
  int packUpDownMaps( buffer_unit_type * & buffer,
                      ElementReferenceAccessor< array1d< localIndex > > const & packList ) const;

  /**
   * @brief Unpack element-to-node and element-to-face maps.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to pack
   * @param overwriteMap flag to indicate whether to overwrite the local map
   * @return the size of data packed.
   */
  int unpackUpDownMaps( buffer_unit_type const * & buffer,
                        ElementReferenceAccessor< localIndex_array > & packList,
                        bool const overwriteMap );

  /**
   * @brief Get the buffer size needed to pack the set of fractured elements and the map toEmbSurfaces.
   * @param packList list of indices to pack
   * @param fractureRegionName name of the fracture region
   * @return the buffer size needed to pack the data
   */
  int packFracturedElementsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList,
                                 string const fractureRegionName ) const;

  /**
   * @brief Pack set of fractured elements and map toEmbSurfaces to a buffer or get the buffer size.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @param fractureRegionName name of the fracture region
   * @return the size of the data packed
   */
  int packFracturedElements( buffer_unit_type * & buffer,
                             ElementViewAccessor< arrayView1d< localIndex > > const & packList,
                             string const fractureRegionName ) const;

  /**
   * @brief Unpack set of fractured elements and map toEmbSurfaces to a buffer or get the buffer size.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @param fractureRegionName name of the fracture region
   * @return the size of the data unpacked
   */
  int unpackFracturedElements( buffer_unit_type const * & buffer,
                               ElementReferenceAccessor< localIndex_array > & packList,
                               string const fractureRegionName );


private:

  /**
   * @brief Pack a list of wrappers or get the buffer size needed to pack.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of the buffer required to pack the wrappers
   */
  template< bool DO_PACKING >
  int packImpl( buffer_unit_type * & buffer,
                ElementViewAccessor< arrayView1d< localIndex > > const & viewAccessor ) const;

  /**
   * @brief Pack a buffer or get the buffer size.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  template< bool DO_PACKING >
  int packGlobalMapsImpl( buffer_unit_type * & buffer,
                          ElementViewAccessor< arrayView1d< localIndex > > const & viewAccessor ) const;

  /**
   * @brief Pack element-to-node and element-to-face maps to a buffer or get the buffer size.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  template< bool DO_PACKING, typename T >
  int
  packUpDownMapsImpl( buffer_unit_type * & buffer,
                      T const & packList ) const;
  /**
   * @brief Unpack element-to-node and element-to-face maps.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to pack
   * @return the size of the data unpacked
   */
  template< typename T >
  int unpackImpl( buffer_unit_type const * & buffer,
                  T & packList );

  /**
   * @brief Pack set of fractured elements and map toEmbSurfaces to a buffer or get the buffer size.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @param fractureRegionName name of the fracture region
   * @return the size of the data packed
   */
  template< bool DO_PACKING >
  int packFracturedElementsImpl( buffer_unit_type * & buffer,
                                 ElementViewAccessor< arrayView1d< localIndex > > const & packList,
                                 string const fractureRegionName ) const;

  /**
   * @brief Copy constructor.
   */
  ElementRegionManager( const ElementRegionManager & );

  /**
   * @brief Copy assignment operator.
   * @return reference to this object
   */
  ElementRegionManager & operator=( const ElementRegionManager & );
};


template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::constructViewAccessor( string const & viewName, string const & neighborName ) const
{
  ElementViewAccessor< LHS > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase const & elemRegion = getRegion( kReg );
    viewAccessor[kReg].resize( elemRegion.numSubRegions() );

    for( typename dataRepository::indexType kSubReg = 0; kSubReg < elemRegion.numSubRegions(); ++kSubReg )
    {
      Group const * group = &elemRegion.getSubRegion( kSubReg );

      if( !neighborName.empty() )
      {
        group = &group->getGroup( ObjectManagerBase::groupKeyStruct::neighborDataString() ).getGroup( neighborName );
      }

      dataRepository::Wrapper< VIEWTYPE > const * const wrapper = group->getWrapperPointer< VIEWTYPE >( viewName );
      if( wrapper )
      {
        viewAccessor[kReg][kSubReg] = wrapper->reference();
      }
    }
  }
  return viewAccessor;
}


template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::
  constructViewAccessor( string const & viewName, string const & neighborName )
{
  ElementViewAccessor< LHS > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase & elemRegion = getRegion( kReg );
    viewAccessor[kReg].resize( elemRegion.numSubRegions() );

    for( typename dataRepository::indexType kSubReg = 0; kSubReg < elemRegion.numSubRegions(); ++kSubReg )
    {
      Group * group = &elemRegion.getSubRegion( kSubReg );

      if( !neighborName.empty() )
      {
        group = &group->getGroup( ObjectManagerBase::groupKeyStruct::neighborDataString() ).getGroup( neighborName );
      }

      dataRepository::Wrapper< VIEWTYPE > * const wrapper = group->getWrapperPointer< VIEWTYPE >( viewName );
      if( wrapper )
      {
        viewAccessor[kReg][kSubReg] = wrapper->reference();
      }
    }
  }
  return viewAccessor;
}

template< typename FIELD_TRAIT >
ElementRegionManager::ElementViewAccessor< traits::ViewTypeConst< typename FIELD_TRAIT::type > >
ElementRegionManager::constructFieldAccessor( string const & neighborName ) const
{
  return constructViewAccessor< typename FIELD_TRAIT::type,
                                traits::ViewTypeConst< typename FIELD_TRAIT::type > >( FIELD_TRAIT::key(), neighborName );
}


template< typename T, int NDIM, typename PERM >
ElementRegionManager::ElementViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
ElementRegionManager::
  constructArrayViewAccessor( string const & name, string const & neighborName ) const
{
  return constructViewAccessor< Array< T, NDIM, PERM >,
                                ArrayView< T const, NDIM, getUSD< PERM > >
                                >( name, neighborName );
}

template< typename VIEWTYPE >
ElementRegionManager::ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
ElementRegionManager::
  constructReferenceAccessor( string const & viewName, string const & neighborName ) const
{
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase const & elemRegion = getRegion( kReg );
    viewAccessor[kReg].resize( elemRegion.numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0; kSubReg<elemRegion.numSubRegions(); ++kSubReg )
    {
      Group const * group = &elemRegion.getSubRegion( kSubReg );

      if( !neighborName.empty() )
      {
        group = &group->getGroup( ObjectManagerBase::groupKeyStruct::neighborDataString() ).getGroup( neighborName );
      }

      if( group->hasWrapper( viewName ) )
      {
        viewAccessor[kReg][kSubReg].set( group->getReference< VIEWTYPE >( viewName ) );
      }
    }
  }
  return viewAccessor;
}

template< typename VIEWTYPE >
ElementRegionManager::ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
ElementRegionManager::constructReferenceAccessor( string const & viewName,
                                                  string const & neighborName )
{
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg = 0; kReg < numRegions(); ++kReg )
  {
    ElementRegionBase & elemRegion = getRegion( kReg );
    viewAccessor[kReg].resize( elemRegion.numSubRegions() );

    for( typename dataRepository::indexType kSubReg = 0; kSubReg < elemRegion.numSubRegions(); ++kSubReg )
    {
      Group * group = &elemRegion.getSubRegion( kSubReg );

      if( !neighborName.empty() )
      {
        group = &group->getGroup( ObjectManagerBase::groupKeyStruct::neighborDataString() ).getGroup( neighborName );
      }

      if( group->hasWrapper( viewName ) )
      {
        viewAccessor[kReg][kSubReg].set( group->getReference< VIEWTYPE >( viewName ) );
      }
    }
  }
  return viewAccessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::MaterialViewAccessor< LHS >
ElementRegionManager::
  constructFullMaterialViewAccessor( string const & viewName,
                                     constitutive::ConstitutiveManager const & cm ) const
{
  MaterialViewAccessor< LHS > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase const & elemRegion = getRegion( kReg );
    accessor[kReg].resize( elemRegion.numSubRegions() );

    for( localIndex kSubReg=0; kSubReg<elemRegion.numSubRegions(); ++kSubReg )
    {
      ElementSubRegionBase const & subRegion = elemRegion.getSubRegion( kSubReg );
      dataRepository::Group const & constitutiveGroup = subRegion.getConstitutiveModels();

      accessor[kReg][kSubReg].resize( cm.numSubGroups() );

      for( localIndex matIndex=0; matIndex<cm.numSubGroups(); ++matIndex )
      {
        string const & constitutiveName = cm.getGroup( matIndex ).getName();
        dataRepository::Group const * const constitutiveRelation = constitutiveGroup.getGroupPointer( constitutiveName );
        if( constitutiveRelation != nullptr )
        {
          dataRepository::Wrapper< VIEWTYPE > const * const wrapper = constitutiveRelation->getWrapperPointer< VIEWTYPE >( viewName );
          if( wrapper )
          {
            accessor[kReg][kSubReg][matIndex] = wrapper->reference();
          }
        }
      }
    }
  }
  return accessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::MaterialViewAccessor< LHS >
ElementRegionManager::
  constructFullMaterialViewAccessor( string const & viewName,
                                     constitutive::ConstitutiveManager const & cm )
{
  MaterialViewAccessor< LHS > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase & elemRegion = getRegion( kReg );
    accessor[kReg].resize( elemRegion.numSubRegions() );

    for( localIndex kSubReg=0; kSubReg<elemRegion.numSubRegions(); ++kSubReg )
    {
      ElementSubRegionBase & subRegion = elemRegion.getSubRegion( kSubReg );
      dataRepository::Group & constitutiveGroup = subRegion.getConstitutiveModels();

      accessor[kReg][kSubReg].resize( cm.numSubGroups() );

      for( localIndex matIndex=0; matIndex<cm.numSubGroups(); ++matIndex )
      {
        string const & constitutiveName = cm.getGroup( matIndex ).getName();
        dataRepository::Group * const constitutiveRelation = constitutiveGroup.getGroupPointer( constitutiveName );
        if( constitutiveRelation != nullptr )
        {
          dataRepository::Wrapper< VIEWTYPE > * const wrapper = constitutiveRelation->getWrapperPointer< VIEWTYPE >( viewName );
          if( wrapper )
          {
            accessor[kReg][kSubReg][matIndex] = wrapper->reference();
          }
        }
      }
    }
  }
  return accessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::constructMaterialViewAccessor( string const & viewName,
                                                     arrayView1d< string const > const & regionNames,
                                                     string const & materialKeyName,
                                                     bool const allowMissingViews ) const
{
  ElementViewAccessor< LHS > accessor;

  // Resize the accessor to all regions and subregions
  accessor.resize( numRegions() );
  for( localIndex kReg = 0; kReg < numRegions(); ++kReg )
  {
    accessor[kReg].resize( getRegion( kReg ).numSubRegions() );
  }

  subGroupMap const & regionMap = getRegions();

  // Loop only over regions named and populate according to given material names
  for( localIndex k = 0; k < regionNames.size(); ++k )
  {
    localIndex const er = regionMap.getIndex( regionNames[k] );
    if( er >=0 )
    {
      GEOS_ERROR_IF_EQ_MSG( er, subGroupMap::KeyIndex::invalid_index, "Region not found: " << regionNames[k] );
      ElementRegionBase const & region = getRegion( er );

      region.forElementSubRegionsIndex( [&]( localIndex const esr,
                                             ElementSubRegionBase const & subRegion )
      {
        string const & materialName = subRegion.getReference< string >( materialKeyName );
        dataRepository::Group const & constitutiveRelation = subRegion.getConstitutiveModel( materialName );

        dataRepository::Wrapper< VIEWTYPE > const * const wrapper = constitutiveRelation.getWrapperPointer< VIEWTYPE >( viewName );
        if( wrapper )
        {
          accessor[er][esr] = wrapper->reference();
        }
        else
        {
          GEOS_ERROR_IF( !allowMissingViews,
                         subRegion.getDataContext() <<
                         ": Material " << constitutiveRelation.getDataContext() <<
                         " does not contain " << viewName );
        }
      } );
    }
  }
  return accessor;
}

template< typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::constructMaterialViewAccessor( string const & viewName,
                                                     arrayView1d< string const > const & regionNames,
                                                     string const & materialKeyName,
                                                     bool const allowMissingViews )
{
  ElementViewAccessor< LHS > accessor;

  // Resize the accessor to all regions and subregions
  accessor.resize( numRegions() );
  for( localIndex kReg = 0; kReg < numRegions(); ++kReg )
  {
    accessor[kReg].resize( getRegion( kReg ).numSubRegions() );
  }

  subGroupMap const & regionMap = getRegions();

  // Loop only over regions named and populate according to given material names
  for( localIndex k = 0; k < regionNames.size(); ++k )
  {
    localIndex const er = regionMap.getIndex( regionNames[k] );
    if( er >=0 )
    {
      GEOS_ERROR_IF_EQ_MSG( er, subGroupMap::KeyIndex::invalid_index, "Region not found: " << regionNames[k] );
      ElementRegionBase & region = getRegion( er );

      region.forElementSubRegionsIndex( [&]( localIndex const esr, ElementSubRegionBase & subRegion )
      {
        string const & materialName = subRegion.getReference< string >( materialKeyName );
        dataRepository::Group const & constitutiveRelation = subRegion.getConstitutiveModel( materialName );

        dataRepository::Wrapper< VIEWTYPE > * const wrapper = constitutiveRelation.getWrapperPointer< VIEWTYPE >( viewName );
        if( wrapper )
        {
          accessor[er][esr] = wrapper->reference();
        }
        else
        {
          GEOS_ERROR_IF( !allowMissingViews, region.getDataContext() << ": Material " << materialName
                                                                     << " does not contain " << viewName );
        }
      } );
    }
  }
  return accessor;
}

template< typename FIELD_TRAIT >
ElementRegionManager::ElementViewAccessor< traits::ViewTypeConst< typename FIELD_TRAIT::type > >
ElementRegionManager::constructMaterialFieldAccessor( arrayView1d< string const > const & regionNames,
                                                      arrayView1d< string const > const & materialNames,
                                                      bool const allowMissingViews ) const
{
  return constructMaterialViewAccessor< typename FIELD_TRAIT::type,
                                        traits::ViewTypeConst< typename FIELD_TRAIT::type > >( FIELD_TRAIT::key(),
                                                                                               regionNames,
                                                                                               materialNames,
                                                                                               allowMissingViews );
}

template< typename MATERIAL_TYPE, typename FIELD_TRAIT >
ElementRegionManager::ElementViewAccessor< traits::ViewTypeConst< typename FIELD_TRAIT::type > >
ElementRegionManager::constructMaterialFieldAccessor( bool const allowMissingViews ) const
{
  GEOS_UNUSED_VAR( allowMissingViews );
  return constructMaterialViewAccessor< MATERIAL_TYPE, typename FIELD_TRAIT::type,
                                        traits::ViewTypeConst< typename FIELD_TRAIT::type > >( FIELD_TRAIT::key() );
}


template< typename T, int NDIM, typename PERM >
ElementRegionManager::ElementViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
ElementRegionManager::
  constructMaterialArrayViewAccessor( string const & viewName,
                                      arrayView1d< string const > const & regionNames,
                                      string const & materialKeyName,
                                      bool const allowMissingViews ) const
{
  return constructMaterialViewAccessor< Array< T, NDIM, PERM >, ArrayView< T const, NDIM, getUSD< PERM > > >( viewName,
                                                                                                              regionNames,
                                                                                                              materialKeyName,
                                                                                                              allowMissingViews );
}

template< typename MATERIALTYPE, typename VIEWTYPE, typename LHS >
ElementRegionManager::ElementViewAccessor< LHS >
ElementRegionManager::constructMaterialViewAccessor( string const & viewName ) const
{
  ElementViewAccessor< LHS > accessor( numRegions() );

  // Resize the accessor to all regions and subregions
  for( localIndex er = 0; er < numRegions(); ++er )
  {
    accessor[er].resize( getRegion( er ).numSubRegions() );
  }

  // Loop only over regions named and populate according to given material names
  for( localIndex er = 0; er < numRegions(); ++er )
  {
    ElementRegionBase const & region = getRegion( er );

    region.forElementSubRegionsIndex( [&]( localIndex const esr,
                                           ElementSubRegionBase const & subRegion )
    {
      dataRepository::Group const & constitutiveGroup = subRegion.getConstitutiveModels();

      string materialName;
      constitutiveGroup.forSubGroups< MATERIALTYPE >( [&]( MATERIALTYPE const & constitutiveRelation )
      {
        materialName = constitutiveRelation.getName();
        if( constitutiveRelation.template hasWrapper( viewName ) )  //NOTE (matteo): I have added this check to allow for the view to be
                                                                    // missing. I am not sure this is the default behaviour we want though.
        {
          accessor[er][esr] = constitutiveRelation.template getReference< VIEWTYPE >( viewName );
        }
      } );
    } );
  }
  return accessor;
}

template< typename MATERIALTYPE, typename T, int NDIM, typename PERM >
ElementRegionManager::ElementViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
ElementRegionManager::constructMaterialArrayViewAccessor( string const & viewName ) const
{
  return constructMaterialViewAccessor< MATERIALTYPE, Array< T, NDIM, PERM >, ArrayView< T const, NDIM, getUSD< PERM > > >( viewName );
}

template< typename CONSTITUTIVE_TYPE >
ElementRegionManager::ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
ElementRegionManager::constructFullConstitutiveAccessor( constitutive::ConstitutiveManager const & cm ) const
{
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase const & elemRegion = getRegion( kReg );
    accessor[kReg].resize( elemRegion.numSubRegions() );

    for( localIndex kSubReg=0; kSubReg<elemRegion.numSubRegions(); ++kSubReg )
    {
      ElementSubRegionBase const & subRegion = elemRegion.getSubRegion( kSubReg );
      dataRepository::Group const & constitutiveGroup = subRegion.getConstitutiveModels();
      accessor[kReg][kSubReg].resize( cm.numSubGroups() );

      for( localIndex matIndex=0; matIndex<cm.numSubGroups(); ++matIndex )
      {
        string const & constitutiveName = cm.getGroup( matIndex ).getName();

        CONSTITUTIVE_TYPE * const
        constitutiveRelation = constitutiveGroup.getGroupPointer< CONSTITUTIVE_TYPE >( constitutiveName );
        if( constitutiveRelation != nullptr )
        {
          accessor[kReg][kSubReg][matIndex] = constitutiveRelation;
        }
      }
    }
  }
  return accessor;
}

template< typename CONSTITUTIVE_TYPE >
ElementRegionManager::ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
ElementRegionManager::constructFullConstitutiveAccessor( constitutive::ConstitutiveManager const & cm )
{
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase & elemRegion = getRegion( kReg );
    accessor[kReg].resize( elemRegion.numSubRegions() );

    for( localIndex kSubReg=0; kSubReg<elemRegion.numSubRegions(); ++kSubReg )
    {
      ElementSubRegionBase & subRegion = elemRegion.getSubRegion( kSubReg );
      dataRepository::Group & constitutiveGroup = subRegion.getConstitutiveModels();
      accessor[kReg][kSubReg].resize( cm.numSubGroups() );

      for( localIndex matIndex=0; matIndex<cm.numSubGroups(); ++matIndex )
      {
        string const & constitutiveName = cm.getGroup( matIndex ).getName();

        CONSTITUTIVE_TYPE * const
        constitutiveRelation = constitutiveGroup.getGroupPointer< CONSTITUTIVE_TYPE >( constitutiveName );
        if( constitutiveRelation != nullptr )
        {
          accessor[kReg][kSubReg][matIndex] = constitutiveRelation;
        }
      }
    }
  }
  return accessor;
}

}
#endif /* GEOS_MESH_ELEMENTREGIONMANAGER_HPP */
