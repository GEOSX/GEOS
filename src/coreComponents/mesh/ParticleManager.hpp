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
 * @file ParticleManager.hpp
 */

#ifndef GEOS_MESH_PARTICLEREGIONMANAGER_HPP
#define GEOS_MESH_PARTICLEREGIONMANAGER_HPP

#include "generators/ParticleBlock.hpp"
#include "generators/ParticleBlockManager.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "ParticleRegion.hpp"
#include "ParticleSubRegion.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "dataRepository/ReferenceWrapper.hpp"

namespace geos
{

class MeshManager;
class DomainPartition;

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
    /// @return particle regions group string key
    static constexpr auto particleRegionsGroup() { return "particleRegionsGroup"; }
  };

  /**
   * @brief The ParticleViewAccessor at the ParticleManager level is an array of array of VIEWTYPE.
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ParticleViewAccessor = array1d< array1d< VIEWTYPE > >;

  /**
   * @brief The ParticleViewAccessor at the ParticleManager level is the
   *   type resulting from ParticleViewAccessor< VIEWTYPE >::toNestedView().
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ParticleView = typename ParticleViewAccessor< VIEWTYPE >::NestedViewType;

  /**
   * @brief The ParticleViewAccessor at the ParticleManager level is the
   *   type resulting from ParticleViewAccessor< VIEWTYPE >::toNestedViewConst().
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ParticleViewConst = typename ParticleViewAccessor< VIEWTYPE >::NestedViewTypeConst;

  /**
   * @brief The ParticleViewAccessor at the ParticleManager level is a 2D array of ReferenceWrapper around VIEWTYPE.
   * @tparam VIEWTYPE data type
   */
  template< typename VIEWTYPE >
  using ParticleReferenceAccessor = array1d< array1d< ReferenceWrapper< VIEWTYPE > > >;

  /**
   * @brief The MaterialViewAccessor at the ParticleManager level is a 3D array of VIEWTYPE.
   * @tparam VIEWTYPE data type
   * var[particleRegionIndex][particleSubRegionIndex][materialIndexInRegion]
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
  static string catalogName()
  { return "ParticleManager"; }

  /**
   * @brief Virtual access to catalogName()
   * @return string that contains the catalog name used to register/lookup this class in the object catalog
   */
  virtual string getCatalogName() const override final
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
   * @brief Get the number of particles in all ParticleSubRegions of type T.
   * @return number of particles
   */
  template< typename T = ParticleSubRegionBase >
  localIndex getNumberOfParticles() const
  {
    localIndex numParticle = 0;
    this->forParticleSubRegions< T >( [&]( ParticleSubRegionBase const & particleSubRegion )
    {
      numParticle += particleSubRegion.size();
    } );
    return numParticle;
  }

//  void Initialize(  ){}

  /**
   * @brief Generate the mesh.
   * @param [in] particleBlockManager pointer to the ParticleBlockManager
   */
  void generateMesh( ParticleBlockManagerABC & particleBlockManager );

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
   * @brief Set the number of particles for a set of particle regions.
   * @param numParticles list of the new particle numbers
   * @param regionNames list of the particle region names
   * @param particleTypes list of the particle types
   */
  void resize( integer_array const & numParticles,
               string_array const & regionNames,
               string_array const & particleTypes );

  /**
   * @brief Set the maximum local and global index.
   */
  virtual void setMaxGlobalIndex() override final;

  /**
   * @brief Get a collection of particle regions
   * @return reference to immutable subGroupMap
   */
  subGroupMap const & getRegions() const
  {
    return this->getGroup( groupKeyStruct::particleRegionsGroup() ).getSubGroups();
  }

  /**
   * @brief Get a collection of particle regions.
   * @return reference to mutable subGroupMap
   */
  subGroupMap & getRegions()
  {
    return this->getGroup( groupKeyStruct::particleRegionsGroup() ).getSubGroups();
  }

  /**
   * @brief Get a particle region.
   * @param key The key of particle region, either name or number.
   * @return Reference to const T.
   */
  template< typename T=ParticleRegionBase, typename KEY_TYPE=void >
  T const & getRegion( KEY_TYPE const & key ) const
  {
    return this->getGroup( groupKeyStruct::particleRegionsGroup() ).getGroup< T >( key );
  }

  /**
   * @brief Get a particle region.
   * @param key The key of the particle region, either name or number.
   * @return Reference to T.
   */
  template< typename T=ParticleRegionBase, typename KEY_TYPE=void >
  T & getRegion( KEY_TYPE const & key )
  {
    return this->getGroup( groupKeyStruct::particleRegionsGroup() ).getGroup< T >( key );
  }

  /**
   * @brief Determines if a ParticleRegion with the input name exists.
   * @tparam T The type of ParticleRegion. May be a specific derived type of ParticleRegionBase.
   * @param name The name/key of the ParticleRegion
   * @return true if the region exists, false if not.
   */
  template< typename T=ParticleRegionBase >
  bool hasRegion( string const & name ) const
  {
    return this->getGroup( groupKeyStruct::particleRegionsGroup() ).hasGroup< T >( name );
  }

  /**
   * @brief Get the number of regions.
   * @return number of regions
   */
  localIndex numRegions() const
  {
    return this->getRegions().size();
  }

  /**
   * @brief Get the number of particle blocks.
   * @return number of particle blocks
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
   * @brief This function is used to launch kernel function over the target particle regions with region type =
   * ParticleRegionBase.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ParticleRegionBase, typename ... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    this->getGroup( groupKeyStruct::particleRegionsGroup() ).forSubGroups< REGIONTYPE, REGIONTYPES... >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the target particle regions with region type =
   * ParticleRegionBase.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
   * @param lambda kernel function
   */
  template< typename REGIONTYPE = ParticleRegionBase, typename ... REGIONTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    this->getGroup( groupKeyStruct::particleRegionsGroup() ).forSubGroups< REGIONTYPE, REGIONTYPES... >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over all the types of particle regions.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forParticleRegionsComplete( LAMBDA lambda ) const
  {
    forParticleRegionsComplete< ParticleRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over all the types of particle regions.
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
   * @brief This const function is used to launch kernel function over the specified target particle regions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda ) const
  {
    forParticleRegionsComplete< ParticleRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target particle regions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA lambda )
  {
    forParticleRegionsComplete< ParticleRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target particle regions with region type =
   * specified particle region types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
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
   * @brief This const function is used to launch kernel function over the specified target particle regions with region
   * type = specified particle region types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
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
   * @brief This function is used to launch kernel function over the particle subregions of all the subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forParticleSubRegions( LAMBDA && lambda )
  {
    forParticleSubRegions< ParticleSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the particle subregions of all the subregion
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
   * @brief This function is used to launch kernel function over the specified target particle subregions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forParticleSubRegions< ParticleSubRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target particle subregions.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleSubRegions( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forParticleSubRegions< ParticleSubRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the particle subregions of the specified subregion
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
   * @brief This const function is used to launch kernel function over the particle subregions of the specified subregion
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
   * @brief This function is used to launch kernel function over the specified target particle subregions with the
   * specified subregion types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
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
   * @brief This const function is used to launch kernel function over the specified target particle subregions with the
   * specified subregion types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
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
   * @brief This const function is used to launch kernel function over the particle subregions of all subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forParticleSubRegionsComplete( LAMBDA && lambda ) const
  {
    forParticleSubRegionsComplete< ParticleSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the particle subregions of all subregion types.
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forParticleSubRegionsComplete( LAMBDA && lambda )
  {
    forParticleSubRegionsComplete< ParticleSubRegion >( std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over the specified target particle subregions
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda )
  {
    forParticleSubRegionsComplete< ParticleSubRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This const function is used to launch kernel function over the specified target particle subregions
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
   * @param lambda kernel function
   */
  template< typename LOOKUP_CONTAINER, typename LAMBDA >
  void forParticleSubRegionsComplete( LOOKUP_CONTAINER const & targetRegions, LAMBDA && lambda ) const
  {
    forParticleSubRegionsComplete< ParticleSubRegion >( targetRegions, std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief This function is used to launch kernel function over all the particle subregions that can be casted to one of
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
   * @brief This const function is used to launch kernel function over all the particle subregions that can be casted to
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
   * @brief This function is used to launch kernel function over the specified target particle subregions that can be
   * casted to one of the specified subregion types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
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
   * @brief This const function is used to launch kernel function over the specified target particle subregions that can
   * be casted to one of the specified subregion types.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetRegions target particle region names or indices
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
   * @brief This is a const function to construct a ParticleViewAccessor to access the data registered on the mesh.
   * @tparam FIELD_TRAIT field type
   * @param neighborName neighbor data name
   * @return ParticleViewAccessor that contains traits::ViewTypeConst< typename TRAIT::type > data
   */
  template< typename FIELD_TRAIT >
  ParticleViewAccessor< traits::ViewTypeConst< typename FIELD_TRAIT::type > >
  constructFieldAccessor( string const & neighborName = string() ) const;

  /**
   * @brief This is a const function to construct a ParticleViewAccessor to access the data registered on the mesh.
   * @tparam VIEWTYPE data type
   * @param name view name of the data
   * @param neighborName neighbor data name
   * @return ParticleViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ParticleViewAccessor< LHS >
  constructViewAccessor( string const & name, string const & neighborName = string() ) const;

  /**
   * @brief This is a function to construct a ParticleViewAccessor to access the data registered on the mesh.
   * @tparam VIEWTYPE data type
   * @param name view name of the data
   * @param neighborName neighbor data name
   * @return ParticleViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ParticleViewAccessor< LHS >
  constructViewAccessor( string const & name, string const & neighborName = string() );

  /**
   * @brief This is a function to construct a ParticleViewAccessor to access array data registered on the mesh.
   * @tparam T data type
   * @tparam NDIM number of array dimensions
   * @tparam PERM layout permutation sequence type
   * @param name view name of the data
   * @param neighborName neighbor data name
   * @return ParticleViewAccessor that contains ArrayView<T const, NDIM> of data
   */
  template< typename T, int NDIM, typename PERM = defaultLayout< NDIM > >
  ParticleViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
  constructArrayViewAccessor( string const & name, string const & neighborName = string() ) const;

  /**
   * @brief This is a const function to construct a ParticleViewAccessor to access the data registered on the mesh.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param neighborName neighbor data name
   * @return ParticleViewAccessor that contains pointers to wrapped VIEWTYPE data
   */
  template< typename VIEWTYPE >
  ParticleViewAccessor< ReferenceWrapper< VIEWTYPE > >
  constructReferenceAccessor( string const & viewName, string const & neighborName = string() ) const;

  /**
   * @brief This is a function to construct a ParticleViewAccessor to access the data registered on the mesh.
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @param neighborName neighbor data name
   * @return ParticleViewAccessor that contains pointers to wrapped VIEWTYPE data
   */
  template< typename VIEWTYPE >
  ParticleViewAccessor< ReferenceWrapper< VIEWTYPE > >
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
   * @return ParticleViewAccessor that contains traits::ViewTypeConst< typename FIELD_TRAIT::type > data
   */
  template< typename FIELD_TRAIT >
  ParticleViewAccessor< traits::ViewTypeConst< typename FIELD_TRAIT::type > >
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
   * @return ParticleViewAccessor that contains traits::ViewTypeConst< typename TRAIT::type > data
   */
  template< typename MATERIAL_TYPE, typename FIELD_TRAIT >
  ParticleViewAccessor< traits::ViewTypeConst< typename FIELD_TRAIT::type > >
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
   * @return ParticleViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ParticleViewAccessor< LHS >
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
   * @return ParticleViewAccessor that contains VIEWTYPE data
   */
  template< typename VIEWTYPE, typename LHS=VIEWTYPE >
  ParticleViewAccessor< LHS >
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
  ParticleViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
  constructMaterialArrayViewAccessor( string const & viewName,
                                      arrayView1d< string const > const & regionNames,
                                      string const & materialKeyName,
                                      bool const allowMissingViews = false ) const;

  /**
   * @brief Construct a const view accessor to material data for specified material type.
   * @tparam MATERIALTYPE base type of material model
   * @tparam VIEWTYPE data type
   * @param viewName view name of the data
   * @return ParticleViewAccessor that contains VIEWTYPE data. Empty views are returned
   *         for subregions that don't contain a model derived from MODELTYPE.
   */
  template< typename MATERIALTYPE, typename VIEWTYPE, typename LHS=VIEWTYPE >
  ParticleViewAccessor< LHS >
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
  ParticleViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
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

  /**
   * @brief Get the buffer size needed to pack a list of wrappers.
   * @param wrapperNames list of wrapper names
   * @param packList list of indices to pack
   * @return the size of the buffer required to pack the wrappers
   */
  int PackSize( string_array const & wrapperNames,
                ParticleViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack a list of wrappers to a buffer.
   * @param buffer pointer to the buffer to be packed
   * @param wrapperNames list of wrapper names
   * @param packList list of indices to pack
   * @return the size of data packed to the buffer
   */
  int Pack( buffer_unit_type * & buffer,
            string_array const & wrapperNames,
            ParticleViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /// @copydoc dataRepository::Group::unpack
  using ObjectManagerBase::unpack;

  /**
   * @brief Unpack a buffer.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to unpack
   * @return the size of data unpacked
   */
  int Unpack( buffer_unit_type const * & buffer,
              ParticleViewAccessor< arrayView1d< localIndex > > & packList );

  /**
   * @brief Unpack a buffer.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to unpack
   * @return the size of data unpacked.
   */
  int Unpack( buffer_unit_type const * & buffer,
              ParticleReferenceAccessor< array1d< localIndex > > & packList );

  /**
   * @brief Get the size of the buffer to be packed.
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  int PackGlobalMapsSize( ParticleViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Pack a buffer.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  int PackGlobalMaps( buffer_unit_type * & buffer,
                      ParticleViewAccessor< arrayView1d< localIndex > > const & packList ) const;

  /**
   * @brief Unpack a buffer.
   * @param buffer pointer to the buffer to be unpacked
   * @param packList list of indices to pack
   * @return the size of the data unpacked
   */
  int UnpackGlobalMaps( buffer_unit_type const * & buffer,
                        ParticleViewAccessor< ReferenceWrapper< localIndex_array > > & packList );

  /**
   * @brief Updates the globalToLocal and localToGlobal maps
   */
  void updateMaps();


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
                   ParticleViewAccessor< arrayView1d< localIndex > > const & viewAccessor ) const;

  /**
   * @brief Pack a buffer or get the buffer size.
   * @param buffer pointer to the buffer to be packed
   * @param packList list of indices to pack
   * @return the size of the data packed
   */
  template< bool DOPACK >
  int PackGlobalMapsPrivate( buffer_unit_type * & buffer,
                             ParticleViewAccessor< arrayView1d< localIndex > > const & viewAccessor ) const;

  /**
   * @brief Unpack particle-to-node and particle-to-face maps.
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

};


template< typename VIEWTYPE, typename LHS >
ParticleManager::ParticleViewAccessor< LHS >
ParticleManager::constructViewAccessor( string const & viewName, string const & neighborName ) const
{
  ParticleViewAccessor< LHS > viewAccessor;
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
ParticleManager::ParticleViewAccessor< LHS >
ParticleManager::constructViewAccessor( string const & viewName, string const & neighborName )
{
  ParticleViewAccessor< LHS > viewAccessor;
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

template< typename FIELD_TRAIT >
ParticleManager::ParticleViewAccessor< traits::ViewTypeConst< typename FIELD_TRAIT::type > >
ParticleManager::constructFieldAccessor( string const & neighborName ) const
{
  return constructViewAccessor< typename FIELD_TRAIT::type,
                                traits::ViewTypeConst< typename FIELD_TRAIT::type > >( FIELD_TRAIT::key(), neighborName );
}

template< typename T, int NDIM, typename PERM >
ParticleManager::ParticleViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
ParticleManager::constructArrayViewAccessor( string const & name, string const & neighborName ) const
{
  return constructViewAccessor< Array< T, NDIM, PERM >,
                                ArrayView< T const, NDIM, getUSD< PERM > >
                                >( name, neighborName );
}

template< typename VIEWTYPE >
ParticleManager::ParticleViewAccessor< ReferenceWrapper< VIEWTYPE > >
ParticleManager::constructReferenceAccessor( string const & viewName, string const & neighborName ) const
{
  ParticleViewAccessor< ReferenceWrapper< VIEWTYPE > > viewAccessor;
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
ParticleManager::ParticleViewAccessor< ReferenceWrapper< VIEWTYPE > >
ParticleManager::constructReferenceAccessor( string const & viewName, string const & neighborName )
{
  ParticleViewAccessor< ReferenceWrapper< VIEWTYPE > > viewAccessor;
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
ParticleManager::constructFullMaterialViewAccessor( string const & viewName,
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
ParticleManager::constructFullMaterialViewAccessor( string const & viewName,
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
ParticleManager::ParticleViewAccessor< LHS >
ParticleManager::constructMaterialViewAccessor( string const & viewName,
                                                arrayView1d< string const > const & regionNames,
                                                string const & materialKeyName,
                                                bool const allowMissingViews ) const
{
  ParticleViewAccessor< LHS > accessor;

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
      ParticleRegionBase const & region = getRegion( er );

      region.forParticleSubRegionsIndex( [&]( localIndex const esr,
                                              ParticleSubRegionBase const & subRegion )
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
          GEOS_ERROR_IF( !allowMissingViews, "Material " << materialKeyName[k] << " does not contain " << viewName );
        }
      } );
    }
  }
  return accessor;
}

template< typename VIEWTYPE, typename LHS >
ParticleManager::ParticleViewAccessor< LHS >
ParticleManager::constructMaterialViewAccessor( string const & viewName,
                                                arrayView1d< string const > const & regionNames,
                                                string const & materialKeyName,
                                                bool const allowMissingViews )
{
  ParticleViewAccessor< LHS > accessor;

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
      ParticleRegionBase & region = getRegion( er );

      region.forParticleSubRegionsIndex( [&]( localIndex const esr, ParticleSubRegionBase & subRegion )
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
          GEOS_ERROR_IF( !allowMissingViews, "Material " << materialName << " does not contain " << viewName );
        }
      } );
    }
  }
  return accessor;
}

template< typename FIELD_TRAIT >
ParticleManager::ParticleViewAccessor< traits::ViewTypeConst< typename FIELD_TRAIT::type > >
ParticleManager::constructMaterialFieldAccessor( arrayView1d< string const > const & regionNames,
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
ParticleManager::ParticleViewAccessor< traits::ViewTypeConst< typename FIELD_TRAIT::type > >
ParticleManager::constructMaterialFieldAccessor( bool const allowMissingViews ) const
{
  GEOS_UNUSED_VAR( allowMissingViews );
  return constructMaterialViewAccessor< MATERIAL_TYPE, typename FIELD_TRAIT::type,
                                        traits::ViewTypeConst< typename FIELD_TRAIT::type > >( FIELD_TRAIT::key() );
}

template< typename T, int NDIM, typename PERM >
ParticleManager::ParticleViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
ParticleManager::constructMaterialArrayViewAccessor( string const & viewName,
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
ParticleManager::ParticleViewAccessor< LHS >
ParticleManager::constructMaterialViewAccessor( string const & viewName ) const
{
  ParticleViewAccessor< LHS > accessor( numRegions() );

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
ParticleManager::ParticleViewAccessor< ArrayView< T const, NDIM, getUSD< PERM > > >
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
#endif /* GEOS_MESH_PARTICLEREGIONMANAGER_HPP */
