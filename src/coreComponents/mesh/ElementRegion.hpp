/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ElementManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 14, 2010
 */

#ifndef ELEMENTREGION_H
#define ELEMENTREGION_H

#include "managers/ObjectManagerBase.hpp"
#include "FaceManager.hpp"
#include "CellBlockSubRegion.hpp"
#include "FaceCellSubRegion.hpp"

namespace geosx
{

//template< typename CELLTYPE, typename ... CELLTYPES, typename LAMBDA >
//constexpr static bool applyLambdaToCellBlocks( ManagedGroup * const cellSubRegion, LAMBDA&& lambda )
//{
//  bool rval = false;
//
//  CELLTYPE * const subRegion = dynamic_cast<CELLTYPE *>( cellSubRegion );
//  if( subRegion!= nullptr )
//  {
//    lambda( subRegion );
//    rval = true;
//  }
//  else
//  {
//    rval = applyLambdaToCellBlocks< CELLTYPES..., LAMBDA >( cellSubRegion, std::forward<LAMBDA>(lambda) );
//  }
//
//  return rval;
//}
//
//template< typename CELLTYPE, typename LAMBDA >
//constexpr static bool applyLambdaToCellBlocks( ManagedGroup * const cellSubRegion, LAMBDA&& lambda )
//{
//  bool rval = false;
//  CELLTYPE * const subRegion = dynamic_cast<CELLTYPE *>( cellSubRegion );
//  if( subRegion!= nullptr )
//  {
//    lambda( subRegion );
//    rval=true;
//  }
//
//  return rval;
//}
//
//template< typename... CELLTYPES, typename LAMBDA >
//void forSomeCellBlocks( LAMBDA && lambda ) const
//{
//  ManagedGroup const * cellBlockSubRegions = this->GetGroup(viewKeyStruct::cellBlockSubRegions);
//
//  for( auto const & subGroupIter : cellBlockSubRegions->GetSubGroups() )
//  {
//    bool isNull =
//    !applyLambdaToCellBlocks<CELLTYPES...>( subGroupIter.second, [&]( auto * const subRegion )
//    {
//      lambda( subRegion );
//    });
//    GEOS_ERROR_IF( isNull, "subRegion "<<subGroupIter.second->getName()<<" is can not be casted to any "
//                   "types specified in parameter pack.");
//  }
//}






class StableTimeStep;


/**
 * Class to manage the data stored at the element level.
 */
class ElementRegion : public ObjectManagerBase
{
public:

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static const string CatalogName()
  { return "ElementRegion"; }

  virtual const string getCatalogName() const override final
  { return ElementRegion::CatalogName(); }


  ///@}


  ElementRegion() = delete;

  ElementRegion( string const & name, ManagedGroup * const parent );


  ElementRegion(const ElementRegion& init);

  virtual ~ElementRegion() override;

  void GenerateMesh( ManagedGroup const * const cellBlocks );

  void GenerateFractureMesh( FaceManager const * const faceManager );

  subGroupMap & GetSubRegions()
  {
    return GetGroup(viewKeyStruct::cellBlockSubRegions)->GetSubGroups();
  }

  subGroupMap const & GetSubRegions() const
  {
    return GetGroup(viewKeyStruct::cellBlockSubRegions)->GetSubGroups();
  }

  template< typename CELLTYPE=CellBase >
  CELLTYPE const * GetSubRegion( string const & regionName ) const
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetGroup<CELLTYPE>(regionName);
  }
  template< typename CELLTYPE=CellBase >
  CELLTYPE * GetSubRegion( string const & regionName )
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetGroup<CELLTYPE>(regionName);
  }

  template< typename CELLTYPE=CellBase >
  CELLTYPE const * GetSubRegion( localIndex const & index ) const
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetGroup<CELLTYPE>(index);
  }
  template< typename CELLTYPE=CellBase >
  CELLTYPE * GetSubRegion( localIndex const & index )
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetGroup<CELLTYPE>(index);
  }

  localIndex numSubRegions() const
  {
    return this->GetGroup(viewKeyStruct::cellBlockSubRegions)->GetSubGroups().size();
  }





  template< typename CELLTYPE1=CellBlockSubRegion, typename CELLTYPE2=CELLTYPE1, typename LAMBDA >
  void forCellBlocks( LAMBDA && lambda ) const
  {
    ManagedGroup const * cellBlockSubRegions = this->GetGroup(viewKeyStruct::cellBlockSubRegions);

    for( auto const & subGroupIter : cellBlockSubRegions->GetSubGroups() )
    {
      bool isNull = true;
      CELLTYPE1 const * const subRegion1 = dynamic_cast<CELLTYPE1 const *>( subGroupIter.second );
      if( subRegion1 != nullptr )
      {
        lambda( subRegion1 );
        isNull = false;
      }
      else
      {
        CELLTYPE2 const * const subRegion2 = dynamic_cast<CELLTYPE2 const *>( subGroupIter.second );
        if( subRegion2 != nullptr )
        {
          lambda( subRegion2 );
          isNull = false;
        }
      }
      GEOS_ERROR_IF( isNull, "subRegion "<<subGroupIter.second->getName()<<" is not of a valid type.");
    }
  }

  template< typename CELLTYPE1=CellBlockSubRegion, typename CELLTYPE2=CELLTYPE1, typename LAMBDA >
  void forCellBlocks( LAMBDA && lambda )
  {
    ManagedGroup * cellBlockSubRegions = this->GetGroup(viewKeyStruct::cellBlockSubRegions);

    for( auto & subGroupIter : cellBlockSubRegions->GetSubGroups() )
    {
      bool isNull = true;
      CELLTYPE1 * const subRegion1 = dynamic_cast<CELLTYPE1 *>( subGroupIter.second );
      if( subRegion1 != nullptr )
      {
        lambda( subRegion1 );
        isNull = false;
      }
      else
      {
        CELLTYPE2 * const subRegion2 = dynamic_cast<CELLTYPE2 *>( subGroupIter.second );
        if( subRegion2 != nullptr )
        {
          lambda( subRegion2 );
          isNull = false;
        }
      }
      GEOS_ERROR_IF( isNull, "subRegion "<<subGroupIter.second->getName()<<" is not of a valid type.");
    }
  }


  template< typename CELLTYPE1=CellBlockSubRegion, typename CELLTYPE2=CELLTYPE1, typename LAMBDA >
  void forCellBlocksIndex( LAMBDA && lambda ) const
  {
    for( localIndex esr=0 ;  esr<this->numSubRegions() ; ++esr )
    {
      CellBase const * const subRegion = this->GetSubRegion(esr);

      bool isNull = true;
      CELLTYPE1 const * const subRegion1 = dynamic_cast<CELLTYPE1  const*>( subRegion );
      if( subRegion1 != nullptr )
      {
        lambda( esr, subRegion1 );
        isNull = false;
      }
      else
      {
        CELLTYPE2 const * const subRegion2 = dynamic_cast<CELLTYPE2 const *>( subRegion );
        if( subRegion2 != nullptr )
        {
          lambda( esr, subRegion2 );
          isNull = false;
        }
      }
      GEOS_ERROR_IF( isNull, "subRegion "<<subRegion->getName()<<" is not of a valid type.");
    }
  }

  template< typename CELLTYPE1=CellBlockSubRegion, typename CELLTYPE2=CELLTYPE1, typename LAMBDA >
  void forCellBlocksIndex( LAMBDA && lambda )
  {
    for( localIndex esr=0 ;  esr<this->numSubRegions() ; ++esr )
    {
      CellBase * const subRegion = this->GetSubRegion(esr);

      bool isNull = true;
      CELLTYPE1 * const subRegion1 = dynamic_cast<CELLTYPE1 *>( subRegion );
      if( subRegion1 != nullptr )
      {
        lambda( esr, subRegion1 );
        isNull = false;
      }
      else
      {
        CELLTYPE2 * const subRegion2 = dynamic_cast<CELLTYPE2 *>( subRegion );
        if( subRegion2 != nullptr )
        {
          lambda( esr, subRegion2 );
          isNull = false;
        }
      }
      GEOS_ERROR_IF( isNull, "subRegion "<<subRegion->getName()<<" is not of a valid type.");
    }
  }


  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto materialListString = "materialList";
    static constexpr auto fractureSetString = "fractureSet";
    static constexpr auto cellBlockSubRegions = "cellBlockSubRegions";
    static constexpr auto cellBlockSubRegionNames = "cellBlocks";

  } m_regionViewKeys;

  string_array & getMaterialList() {return m_materialList;}
  string_array const & getMaterialList() const {return m_materialList;}

protected:
  virtual void PostProcessInput() override;

private:

  ElementRegion& operator=(const ElementRegion& rhs);

  string_array m_cellBlockNames;
  string_array m_fractureSetNames;
  string_array m_materialList;
  string m_numericalMethod;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



}



#endif /* ELEMENTOBJECTT_H_ */
