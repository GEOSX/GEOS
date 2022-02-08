/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ParticleManager.hpp
 */

#ifndef GEOSX_MESH_PARTICLEREGIONMANAGER_HPP
#define GEOSX_MESH_PARTICLEREGIONMANAGER_HPP

#include "ParticleBlock.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "ParticleRegion.hpp"
#include "ParticleSubRegion.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "dataRepository/ReferenceWrapper.hpp"

namespace geosx
{

class MeshManager;

/**
 * @class ParticleManager
 * @brief The ParticleManager class provides an interface to ObjectManagerBase in order to manage ParticleRegion
 * data
 */
class ParticleManager : public ObjectManagerBase
{
public:

  /**
   * @brief Group key associated with particleRegionsGroup.
   */
  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    /// @return element regions group string key
    static constexpr auto particleRegionsGroup() { return "particleRegionsGroup"; }
  };

  /**
   * Limit on max number of nodes for each element
   */
  constexpr static int maxNumNodesPerElem = 8;

  /**
   * @brief The ElementViewAccessor at the ParticleManager level is an array of array of VIEWTYPE.
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ElementViewAccessor = array1d< array1d< VIEWTYPE > >;

  /**
   * @brief The ElementViewAccessor at the ParticleManager level is the
   *   type resulting from ElementViewAccessor< VIEWTYPE >::toNestedView().
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ElementView = typename ElementViewAccessor< VIEWTYPE >::NestedViewType;

  /**
   * @brief The ElementViewAccessor at the ParticleManager level is the
   *   type resulting from ElementViewAccessor< VIEWTYPE >::toNestedViewConst().
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ElementViewConst = typename ElementViewAccessor< VIEWTYPE >::NestedViewTypeConst;

  /**
   * @brief The ElementViewAccessor at the ParticleManager level is a 2D array of ReferenceWrapper around VIEWTYPE.
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ElementReferenceAccessor = array1d< array1d< ReferenceWrapper< VIEWTYPE > > >;

  /**
   * @brief The MaterialViewAccessor at the ParticleManager level is a 3D array of VIEWTYPE.
   * @tparam VIEWTYPE data type
   * var[particleRegionIndex][elementSubRegionIndex][materialIndexInRegion]
   */
  template< typename VIEWTYPE >
  using MaterialViewAccessor = array1d< array1d< array1d< VIEWTYPE > > >;

  /**
   * @brief The ConstitutiveRelationAccessor at the ParticleManager level is a 3D array of CONSTITUTIVE_TYPE
   * @tparam CONSTITUTIVE_TYPE constitutive type
   */
  template< typename CONSTITUTIVE_TYPE >
  using ConstitutiveRelationAccessor = array1d< array1d< array1d< CONSTITUTIVE_TYPE * > > >;

  /**
   * @brief The function is to return the name of the ParticleManager in the object catalog
   * @return string that contains the catalog name used to register/lookup this class in  the object catalog
   */
  static const string catalogName()
  { return "ZoneManager"; }

  /**
   * @brief Virtual access to catalogName()
   * @return string that contains the catalog name used to register/lookup this class in the object catalog
   */
  virtual const string getCatalogName() const override final
  { return ParticleManager::catalogName(); }

  /**
   * @brief Constructor.
   * @param [in] name the name of this ObjectManager
   * @param [in] parent the parent Group
   */
  ParticleManager( string const & name, Group * const parent );

  /**
   * @brief Destructor
   */
  virtual ~ParticleManager() override;

  /**
   * @brief Get the number of elements in all ParticleSubRegions of type T.
   * @return number of elements
   */
  template< typename T = ParticleSubRegionBase >
  localIndex getNumberOfElements() const
  {
    localIndex numElem = 0;
    this->forParticleSubRegions< T >( [&]( ParticleSubRegionBase const & particleBlock )
    {
      numElem += particleBlock.size();
    } );
    return numElem;
  }

//  void Initialize(  ){}

  /**
   * @brief Generate the mesh.
   * @param [in] particleBlockManager pointer to the ParticleBlockManager
   */
  void generateMesh( Group & particleBlockManager );

  /**
   * @brief Generate the aggregates.
   * @param [in] faceManager pointer to the FaceManager
   * @param [in] nodeManager pointer to the NodeManager
   */
//  void generateAggregates( FaceManager const & faceManager, NodeManager const & nodeManager );

  /**
   * @brief Create a new ParticleRegion object as a child of this group.
   * @param childKey catalog key of the new ParticleRegion derived type to create
   * @param childName name of the new ParticleRegion object
   * @return pointer to the created ParticleRegion object
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
   * @param numParticles list of the new element numbers
   * @param regionNames list of the element region names
   * @param particleTypes list of the element types
   */
  void resize( integer_array const & numParticles,
               string_array const & regionNames,
               string_array const & particleTypes );

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
    return this->getGroup( groupKeyStruct::particleRegionsGroup() ).getSubGroups();
  }

  /**
   * @brief Get a collection of element regions.
   * @return reference to mutable subGroupMap
   */
  subGroupMap & getRegions()
  {
    return this->getGroup( groupKeyStruct::particleRegionsGroup() ).getSubGroups();
  }

  /**
   * @brief Get a element region.
   * @param key The key of element region, either name or number.
   * @return Reference to const T.
   */
  template< typename T=ParticleRegionBase, typename KEY_TYPE=void >
  T const & getRegion( KEY_TYPE const & key ) const
  {
    return this->getGroup( groupKeyStruct::particleRegionsGroup() ).getGroup< T >( key );
  }

  /**
   * @brief Get a element region.
   * @param key The key of the element region, either name or number.
   * @return Reference to T.
   */
  template< typename T=ParticleRegionBase, typename KEY_TYPE=void >
  T & getRegion( KEY_TYPE const & key )
  {
    return this->getGroup( groupKeyStruct::particleRegionsGroup() ).getGroup< T >( key );
  }

  template< typename T=ParticleRegionBase >
  bool hasRegion( string const & name ) const
  {
    return this->getGroup( groupKeyStruct::particleRegionsGroup() ).hasGroup< T >( name );
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
   * @brief Get number of the cell blocks.
   * @return number of the cell blocks
   */
  localIndex numParticleBlocks() const;

  /**
   * @brief This function is used to launch kernel function over all the particle regions with region type =
   * ParticleRegionBase.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ParticleRegionBase, typename ... REGIONTYPES, typename LAMBDA >
  void forParticleRegions( LAMBDA && lambda )
  {
    this->getGroup( groupKeyStruct::particleRegionsGroup() ).forSubGroups< REGIONTYPE, REGIONTYPES... >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over all the particle regions with region type =
   * ParticleRegionBase.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ParticleRegionBase, typename ... REGIONTYPES, typename LAMBDA >
  void forParticleRegions( LAMBDA && lambda ) const
  {
    this->getGroup( groupKeyStruct::particleRegionsGroup() ).forSubGroups< REGIONTYPE, REGIONTYPES... >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the target element regions with region type =
   * ParticleRegionBase.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ParticleRegionBase, typename ... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    this->getGroup( groupKeyStruct::particleRegionsGroup() ).forSubGroups< REGIONTYPE, REGIONTYPES... >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the target element regions with region type =
   * ParticleRegionBase.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ParticleRegionBase, typename ... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    this->getGroup( groupKeyStruct::particleRegionsGroup() ).forSubGroups< REGIONTYPE, REGIONTYPES... >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over all the types of element regions.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forParticleRegionsComplete( LAMBDA lambda ) const
  {
    forParticleRegionsComplete< ParticleRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over all the types of element regions.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forParticleRegionsComplete( LAMBDA lambda )
  {
    forParticleRegionsComplete< ParticleRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over all the particle regions that can be casted to one of
   * the specified region types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename REGIONTYPE, typename ... REGIONTYPES, typename LAMBDA >
  void forParticleRegionsComplete( LAMBDA lambda )
  {
    for( localIndex er=0; er<this->numRegions(); ++er )
    {
      ParticleRegionBase & particleRegion = this->getRegion( er );

      Group::applyLambdaToContainer< REGIONTYPE, REGIONTYPES... >( particleRegion, [&]( auto & castedRegion )
      {
        lambda( er, castedRegion );
      } );
    }
  }

  /**
   * @brief This const function is used to launch kernel function over all the particle regions that can be casted to one
   * of the specified region types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename REGIONTYPE, typename ... REGIONTYPES, typename LAMBDA >
  void forParticleRegionsComplete( LAMBDA lambda ) const
  {
    for( localIndex er=0; er<this->numRegions(); ++er )
    {
      ParticleRegionBase const & particleRegion = this->getRegion( er );

      Group::applyLambdaToContainer< REGIONTYPE, REGIONTYPES... >( particleRegion, [&]( auto const & castedRegion )
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
  void forParticleRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda ) const
  {
    forParticleRegionsComplete< ParticleRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element regions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda )
  {
    forParticleRegionsComplete< ParticleRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
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
  void forParticleRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda )
  {
    forParticleRegions< REGIONTYPE, REGIONTYPES... >( targetRegions, [&] ( localIndex const targetIndex,
                                                                          auto & particleRegion )
    {
      lambda( targetIndex, particleRegion.getIndexInParent(), particleRegion );
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
  void forParticleRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda ) const
  {
    forParticleRegions< REGIONTYPE, REGIONTYPES... >( targetRegions, [&] ( localIndex const targetIndex,
                                                                          auto const & particleRegion )
    {
      lambda( targetIndex, particleRegion.getIndexInParent(), particleRegion );
    } );
  }

  /**
   * @brief This function is used to launch kernel function over the element subregions of all the subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forParticleSubRegions( LAMBDA && lambda )
  {
    forParticleSubRegions< ParticleSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the element subregions of all the subregion
   * types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forParticleSubRegions( LAMBDA && lambda ) const
  {
    forParticleSubRegions< ParticleSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element subregions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forParticleSubRegions< ParticleSubRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element subregions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forParticleSubRegions< ParticleSubRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the element subregions of the specified subregion
   * types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forParticleSubRegions( LAMBDA && lambda )
  {
    forParticleSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >(
      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const,
                                                   localIndex const,
                                                   ParticleRegionBase &,
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
  void forParticleSubRegions( LAMBDA && lambda ) const
  {
    forParticleSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >(
      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const,
                                                   localIndex const,
                                                   ParticleRegionBase const &,
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
  void forParticleSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forParticleSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >( targetRegions,
                                                                      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const targetIndex,
                                                                                                                   localIndex const,
                                                                                                                   localIndex const,
                                                                                                                   ParticleRegionBase &,
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
  void forParticleSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forParticleSubRegionsComplete< SUBREGIONTYPE, SUBREGIONTYPES... >( targetRegions,
                                                                      [lambda = std::forward< LAMBDA >( lambda )]( localIndex const targetIndex,
                                                                                                                   localIndex const,
                                                                                                                   localIndex const,
                                                                                                                   ParticleRegionBase const &,
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
  void forParticleSubRegionsComplete( LAMBDA && lambda ) const
  {
    forParticleSubRegionsComplete< ParticleSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the element subregions of all subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forParticleSubRegionsComplete( LAMBDA && lambda )
  {
    forParticleSubRegionsComplete< ParticleSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target element subregions
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forParticleSubRegionsComplete< ParticleSubRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target element subregions
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target element region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forParticleSubRegionsComplete< ParticleSubRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over all the element subregions that can be casted to one of
   * the specified subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename SUBREGIONTYPE, typename ... SUBREGIONTYPES, typename LAMBDA >
  void forParticleSubRegionsComplete( LAMBDA && lambda )
  {
    for( localIndex er=0; er<this->numRegions(); ++er )
    {
      ParticleRegionBase & particleRegion = this->getRegion( er );

      for( localIndex esr=0; esr<particleRegion.numSubRegions(); ++esr )
      {
        ParticleSubRegionBase & subRegion = particleRegion.getSubRegion( esr );

        Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto & castedSubRegion )
        {
          lambda( er, esr, particleRegion, castedSubRegion );
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
  void forParticleSubRegionsComplete( LAMBDA && lambda ) const
  {
    for( localIndex er=0; er<this->numRegions(); ++er )
    {
      ParticleRegionBase const & particleRegion = this->getRegion( er );

      for( localIndex esr=0; esr<particleRegion.numSubRegions(); ++esr )
      {
        ParticleSubRegionBase const & subRegion = particleRegion.getSubRegion( esr );

        Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto const & castedSubRegion )
        {
          lambda( er, esr, particleRegion, castedSubRegion );
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
  void forParticleSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forParticleRegions( targetRegions, [&] ( localIndex const targetIndex, ParticleRegionBase & particleRegion )
    {
      localIndex const er = particleRegion.getIndexInParent();

      if( er>-1 )
      {
        for( localIndex esr=0; esr<particleRegion.numSubRegions(); ++esr )
        {
          ParticleSubRegionBase & subRegion = particleRegion.getSubRegion( esr );

          Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto & castedSubRegion )
          {
            lambda( targetIndex, er, esr, particleRegion, castedSubRegion );
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
  void forParticleSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forParticleRegions( targetRegions, [&] ( localIndex const targetIndex, ParticleRegionBase const & particleRegion )
    {
      localIndex const er = particleRegion.getIndexInParent();

      if( er>-1 )
      {
        for( localIndex esr=0; esr<particleRegion.numSubRegions(); ++esr )
        {
          ParticleSubRegionBase const & subRegion = particleRegion.getSubRegion( esr );

          Group::applyLambdaToContainer< SUBREGIONTYPE, SUBREGIONTYPES... >( subRegion, [&]( auto const & castedSubRegion )
          {
            lambda( targetIndex, er, esr, particleRegion, castedSubRegion );
          } );
        }
      }
    } );
  }


  /**
   * @brief This is a const function to construct a ElementViewAccessor to access the data registered on the mesh.
   * @tparam TRAIT data type
   * @param neighborName neighbor data name
   * @return ElementViewAccessor that contains traits::ViewTypeConst< typename TRAIT::type > data
   */
  template< typename TRAIT >
  ElementViewAccessor< traits::ViewTypeConst< typename TRAIT::type > >
  constructExtrinsicAccessor( string const & neighborName = string() ) const;

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
   * @tparam TRAIT mesh data trait
   * @param regionNames list of region names
   * @param materialNames list of corresponding material names
   * @param allowMissingViews flag to indicate whether it is allowed to miss the specified material data in material
   * list
   * @return ElementViewAccessor that contains traits::ViewTypeConst< typename TRAIT::type > data
   */
  template< typename TRAIT >
  ElementViewAccessor< traits::ViewTypeConst< typename TRAIT::type > >
  constructMaterialExtrinsicAccessor( arrayView1d< string const > const & regionNames,
                                      arrayView1d< string const > const & materialNames,
                                      bool const allowMissingViews = false ) const;

  /**
   * @brief This is a const function to construct a MaterialViewAccessor to access the material data for specified
   * material type.
   * @tparam MATERIALTYPE base type of material model
   * @tparam TRAIT mesh data trait
   * @param allowMissingViews flag to indicate whether it is allowed to miss the specified material data in material
   * list
   * @return ElementViewAccessor that contains traits::ViewTypeConst< typename TRAIT::type > data
   */
  template< typename MATERIALTYPE, typename TRAIT >
  ElementViewAccessor< traits::ViewTypeConst< typename TRAIT::type > >
  constructMaterialExtrinsicAccessor( bool const allowMissingViews = false ) const;


  /**
   * @brief This is a const function to construct a MaterialViewAccessor to access the material data for specified
   * regions/materials.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param regionNames list of region names
   * @param materialNames list of corresponding material names
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
   * @param materialNames list of corresponding material names
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
   * @param materialNames list of corresponding material names
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
   * @brief Get the buffer size needed to pack a list of wrappers.
   * @param wrapperNames list of wrapper names
   * @param packList list of indices to pack
   * @return the size of the buffer required to pack the wrappers
   */
  int PackSize( string_array const & wrapperNames,
                ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack a list of wrappers to a buffer.
   * @param buffer pointer to the buffer to be packed
   * @param wrapperNames list of wrapper names
   * @param packList list of indices to pack
   * @return the size of data packed to the buffer
   */
  int Pack( buffer_unit_type * & buffer,
            string_array const & wrapperNames,
            ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /// @copydoc dataRepository::Group::unpack
  using ObjectManagerBase::unpack;

  /**
   * @brief Unpack a buffer.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to unpack
   * @return the size of data unpacked
   */
  int Unpack( buffer_unit_type const * & buffer,
              ElementViewAccessor< arrayView1d< localIndex > > & packList );

  /**
   * @brief Unpack a buffer.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to unpack
   * @return the size of data unpacked.
   */
  int Unpack( buffer_unit_type const * & buffer,
              ElementReferenceAccessor< array1d< localIndex > > & packList );

  /**
   * @brief Get the size of the buffer to be packed.
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  int PackGlobalMapsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack a buffer.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  int PackGlobalMaps( buffer_unit_type * & buffer,
                      ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Unpack a buffer.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to pack
   * @return the size of the data unpacked
   */
  int UnpackGlobalMaps( buffer_unit_type const * & buffer,
                        ElementViewAccessor< ReferenceWrapper< localIndex_array > > & packList );

  /**
   * @brief Get the buffer size needed to pack element-to-node and element-to-face maps.
   * @param packList list of indices to pack
   * @return the size of data packed.
   */
  int PackUpDownMapsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Get the buffer size needed to pack element-to-node and element-to-face maps.
   * @param packList list of indices to pack
   * @return the size of data packed.
   */
  int PackUpDownMapsSize( ElementReferenceAccessor< array1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack element-to-node and element-to-face maps.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of data packed.
   */
  int PackUpDownMaps( buffer_unit_type * & buffer,
                      ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack element-to-node and element-to-face maps.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of data packed.
   */
  int PackUpDownMaps( buffer_unit_type * & buffer,
                      ElementReferenceAccessor< array1d< localIndex > > const & packList ) const;

  /**
   * @brief Unpack element-to-node and element-to-face maps.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to pack
   * @param overwriteMap flag to indicate whether to overwrite the local map
   * @return the size of data packed.
   */
  int UnpackUpDownMaps( buffer_unit_type const * & buffer,
                        ElementReferenceAccessor< localIndex_array > & packList,
                        bool const overwriteMap );

  /**
   * @name viewKeyStruct/groupKeyStruct
   */
  ///@{

  /**
   *  @brief contains the added view access keys to be bound with class data member.
   *  @struct viewKeyStruct
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// @return String to access the reference position
    static constexpr char const * referencePositionString() { return "ReferencePosition"; }

    /// Accessor to reference position
    dataRepository::ViewKey referencePosition       = { referencePositionString() };
  }
  /// viewKeys
  viewKeys;

  //START_SPHINX_REFPOS_ACCESS
  /**
   * @brief Get the mutable reference position array. This table will contain all the node coordinates.
   * @return reference position array
   */
  array2d< real64, nodes::REFERENCE_POSITION_PERM > & referencePosition() { return m_referencePosition; }

  /**
   * @brief Provide an immutable arrayView of the reference position. This table will contain all the node coordinates.
   * @return an immutable arrayView of the reference position.
   */

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > referencePosition() const
  { return m_referencePosition; }
  //END_SPHINX_REFPOS_ACCESS

private:

  /**
   * @brief Pack a list of wrappers or get the buffer size needed to pack.
   * @param buffer pointer to the buffer to be packed
   * @param wrapperNames list of wrapper names
   * @param packList list of indices to pack
   * @return the size of the buffer required to pack the wrappers
   */
  template< bool DOPACK >
  int PackPrivate( buffer_unit_type * & buffer,
                   string_array const & wrapperNames,
                   ElementViewAccessor< arrayView1d< localIndex > > const & viewAccessor ) const;

  /**
   * @brief Pack a buffer or get the buffer size.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  template< bool DOPACK >
  int PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                             ElementViewAccessor< arrayView1d< localIndex > > const & viewAccessor ) const;

  /**
   * @brief Pack element-to-node and element-to-face maps to a buffer or get the buffer size.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  template< bool DOPACK, typename T >
  int
  packUpDownMapsPrivate( buffer_unit_type * & buffer,
                         T const & packList ) const;
  /**
   * @brief Unpack element-to-node and element-to-face maps.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to pack
   * @return the size of the data unpacked
   */
  template< typename T >
  int unpackPrivate( buffer_unit_type const * & buffer,
                     T & packList );

  /**
   * @brief Copy constructor.
   */
  ParticleManager( const ParticleManager & );

  /**
   * @brief Copy assignment operator.
   * @return reference to this object
   */
  ParticleManager & operator=( const ParticleManager & );

  //START_SPHINX_REFPOS
  /// reference position of the nodes
  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_referencePosition;
  //END_SPHINX_REFPOS

};


template< typename VIEWTYPE, typename LHS >
ParticleManager::ElementViewAccessor< LHS >
ParticleManager::constructViewAccessor( string const & viewName, string const & neighborName ) const
{
  ElementViewAccessor< LHS > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ParticleRegionBase const & particleRegion = getRegion( kReg );
    viewAccessor[kReg].resize( particleRegion.numSubRegions() );

    for( typename dataRepository::indexType kSubReg = 0; kSubReg < particleRegion.numSubRegions(); ++kSubReg )
    {
      Group const * group = &particleRegion.getSubRegion( kSubReg );

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
ParticleManager::ElementViewAccessor< LHS >
ParticleManager::
  constructViewAccessor( string const & viewName, string const & neighborName )
{
  ElementViewAccessor< LHS > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ParticleRegionBase & particleRegion = getRegion( kReg );
    viewAccessor[kReg].resize( particleRegion.numSubRegions() );

    for( typename dataRepository::indexType kSubReg = 0; kSubReg < particleRegion.numSubRegions(); ++kSubReg )
    {
      Group * group = &particleRegion.getSubRegion( kSubReg );

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

template< typename TRAIT >
ParticleManager::ElementViewAccessor< traits::ViewTypeConst< typename TRAIT::type > >
ParticleManager::
  constructExtrinsicAccessor( string const & neighborName ) const
{
  return constructViewAccessor< typename TRAIT::type,
                                traits::ViewTypeConst< typename TRAIT::type > >( TRAIT::key(), neighborName );
}


template< typename T, int NDIM, typename PERM >
ParticleManager::ElementViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
ParticleManager::
  constructArrayViewAccessor( string const & name, string const & neighborName ) const
{
  return constructViewAccessor< Array< T, NDIM, PERM >,
                                ArrayView< T const, NDIM, getUSD< PERM > >
                                >( name, neighborName );
}

template< typename VIEWTYPE >
ParticleManager::ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
ParticleManager::
  constructReferenceAccessor( string const & viewName, string const & neighborName ) const
{
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ParticleRegionBase const & particleRegion = getRegion( kReg );
    viewAccessor[kReg].resize( particleRegion.numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0; kSubReg<particleRegion.numSubRegions(); ++kSubReg )
    {
      Group const * group = &particleRegion.getSubRegion( kSubReg );

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
ParticleManager::ElementViewAccessor< ReferenceWrapper< VIEWTYPE > >
ParticleManager::
  constructReferenceAccessor( string const & viewName, string const & neighborName )
{
  ElementViewAccessor< ReferenceWrapper< VIEWTYPE > > viewAccessor;
  viewAccessor.resize( numRegions() );
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ParticleRegionBase & particleRegion = getRegion( kReg );
    viewAccessor[kReg].resize( particleRegion.numSubRegions() );

    for( typename dataRepository::indexType kSubReg=0; kSubReg<particleRegion.numSubRegions(); ++kSubReg )
    {
      Group * group = &particleRegion.getSubRegion( kSubReg );

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
ParticleManager::MaterialViewAccessor< LHS >
ParticleManager::
  constructFullMaterialViewAccessor( string const & viewName,
                                     constitutive::ConstitutiveManager const & cm ) const
{
  MaterialViewAccessor< LHS > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0; kReg<numRegions(); ++kReg )
  {
    ParticleRegionBase const & particleRegion = getRegion( kReg );
    accessor[kReg].resize( particleRegion.numSubRegions() );

    for( localIndex kSubReg=0; kSubReg<particleRegion.numSubRegions(); ++kSubReg )
    {
      ParticleSubRegionBase const & subRegion = particleRegion.getSubRegion( kSubReg );
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
ParticleManager::MaterialViewAccessor< LHS >
ParticleManager::
  constructFullMaterialViewAccessor( string const & viewName,
                                     constitutive::ConstitutiveManager const & cm )
{
  MaterialViewAccessor< LHS > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0; kReg<numRegions(); ++kReg )
  {
    ParticleRegionBase & particleRegion = getRegion( kReg );
    accessor[kReg].resize( particleRegion.numSubRegions() );

    for( localIndex kSubReg=0; kSubReg<particleRegion.numSubRegions(); ++kSubReg )
    {
      ParticleSubRegionBase & subRegion = particleRegion.getSubRegion( kSubReg );
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
ParticleManager::ElementViewAccessor< LHS >
ParticleManager::constructMaterialViewAccessor( string const & viewName,
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
      GEOSX_ERROR_IF_EQ_MSG( er, subGroupMap::KeyIndex::invalid_index, "Region not found: " << regionNames[k] );
      ParticleRegionBase const & region = getRegion( er );

      region.forParticleSubRegionsIndex( [&]( localIndex const esr,
                                             ParticleSubRegionBase const & subRegion )
      {
        string const & materialName = subRegion.getReference<string>( materialKeyName );
        dataRepository::Group const & constitutiveRelation = subRegion.getConstitutiveModel(materialName);

        dataRepository::Wrapper< VIEWTYPE > const * const wrapper = constitutiveRelation.getWrapperPointer< VIEWTYPE >( viewName );
        if( wrapper )
        {
          accessor[er][esr] = wrapper->reference();
        }
        else
        {
          GEOSX_ERROR_IF( !allowMissingViews, "Material " << materialKeyName[k] << " does not contain " << viewName );
        }
      } );
    }
  }
  return accessor;
}

template< typename VIEWTYPE, typename LHS >
ParticleManager::ElementViewAccessor< LHS >
ParticleManager::constructMaterialViewAccessor( string const & viewName,
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
      GEOSX_ERROR_IF_EQ_MSG( er, subGroupMap::KeyIndex::invalid_index, "Region not found: " << regionNames[k] );
      ParticleRegionBase & region = getRegion( er );

      region.forParticleSubRegionsIndex( [&]( localIndex const esr, ParticleSubRegionBase & subRegion )
      {
        string const & materialName = subRegion.getReference<string>( materialKeyName );
        dataRepository::Group const & constitutiveRelation = subRegion.getConstitutiveModel(materialName);

        dataRepository::Wrapper< VIEWTYPE > * const wrapper = constitutiveRelation.getWrapperPointer< VIEWTYPE >( viewName );
        if( wrapper )
        {
          accessor[er][esr] = wrapper->reference();
        }
        else
        {
          GEOSX_ERROR_IF( !allowMissingViews, "Material " << materialName << " does not contain " << viewName );
        }
      } );
    }
  }
  return accessor;
}

template< typename TRAIT >
ParticleManager::ElementViewAccessor< traits::ViewTypeConst< typename TRAIT::type > >
ParticleManager::
  constructMaterialExtrinsicAccessor( arrayView1d< string const > const & regionNames,
                                      arrayView1d< string const > const & materialNames,
                                      bool const allowMissingViews ) const
{
  return constructMaterialViewAccessor< typename TRAIT::type,
                                        traits::ViewTypeConst< typename TRAIT::type > >( TRAIT::key(),
                                                                                         regionNames,
                                                                                         materialNames,
                                                                                         allowMissingViews );
}

template< typename MATERIALTYPE, typename TRAIT >
ParticleManager::ElementViewAccessor< traits::ViewTypeConst< typename TRAIT::type > >
ParticleManager::
  constructMaterialExtrinsicAccessor( bool const allowMissingViews ) const
{
  GEOSX_UNUSED_VAR(allowMissingViews);
  return constructMaterialViewAccessor< MATERIALTYPE, typename TRAIT::type,
                                        traits::ViewTypeConst< typename TRAIT::type > >( TRAIT::key() );
}


template< typename T, int NDIM, typename PERM >
ParticleManager::ElementViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
ParticleManager::
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
ParticleManager::ElementViewAccessor< LHS >
ParticleManager::constructMaterialViewAccessor( string const & viewName ) const
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
    ParticleRegionBase const & region = getRegion( er );

    region.forParticleSubRegionsIndex( [&]( localIndex const esr,
                                           ParticleSubRegionBase const & subRegion )
    {
      dataRepository::Group const & constitutiveGroup = subRegion.getConstitutiveModels();

      string materialName;
      constitutiveGroup.forSubGroups< MATERIALTYPE >( [&]( MATERIALTYPE const & constitutiveRelation )
      {
        materialName = constitutiveRelation.getName();
        if ( constitutiveRelation.template hasWrapper( viewName ) ) //NOTE (matteo): I have added this check to allow for the view to be missing. I am not sure this is the default behaviour we want though.
        {
          accessor[er][esr] = constitutiveRelation.template getReference< VIEWTYPE >( viewName );
        }
      } );
    } );
  }
  return accessor;
}

template< typename MATERIALTYPE, typename T, int NDIM, typename PERM >
ParticleManager::ElementViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
ParticleManager::constructMaterialArrayViewAccessor( string const & viewName ) const
{
  return constructMaterialViewAccessor< MATERIALTYPE, Array< T, NDIM, PERM >, ArrayView< T const, NDIM, getUSD< PERM > > >( viewName );
}

template< typename CONSTITUTIVE_TYPE >
ParticleManager::ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
ParticleManager::constructFullConstitutiveAccessor( constitutive::ConstitutiveManager const & cm ) const
{
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0; kReg<numRegions(); ++kReg )
  {
    ParticleRegionBase const & particleRegion = getRegion( kReg );
    accessor[kReg].resize( particleRegion.numSubRegions() );

    for( localIndex kSubReg=0; kSubReg<particleRegion.numSubRegions(); ++kSubReg )
    {
      ParticleSubRegionBase const & subRegion = particleRegion.getSubRegion( kSubReg );
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
ParticleManager::ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE >
ParticleManager::constructFullConstitutiveAccessor( constitutive::ConstitutiveManager const & cm )
{
  ConstitutiveRelationAccessor< CONSTITUTIVE_TYPE > accessor;
  accessor.resize( numRegions() );
  for( localIndex kReg=0; kReg<numRegions(); ++kReg )
  {
    ParticleRegionBase & particleRegion = getRegion( kReg );
    accessor[kReg].resize( particleRegion.numSubRegions() );

    for( localIndex kSubReg=0; kSubReg<particleRegion.numSubRegions(); ++kSubReg )
    {
      ParticleSubRegionBase & subRegion = particleRegion.getSubRegion( kSubReg );
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
#endif /* GEOSX_MESH_PARTICLEREGIONMANAGER_HPP */
