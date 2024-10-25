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
 * @file Utilities.cpp
 */

#include "common/MpiWrapper.hpp"
#include "Utilities.hpp"
#include "Group.hpp"


#include <unordered_set>
#include <unordered_map>

namespace geos
{
namespace dataRepository
{

void printMemoryAllocation( Group const & group, integer const indent, real64 const threshold )
{
  // static flag to keep track of whether or not a subgroup is the last one, and thus
  // visually will terminate the tree.
  static bool terminateBranch[64]{};

  // keep track of the address of the groups that have been printed. Not the best
  // way to do this, but doens't require infrastructure changes.
  static std::unordered_set< Group const * > groupsPrinted;
  // initialize the static variables at the head of the tree.
  if( indent==0 )
  {
    terminateBranch[0] = true;
    for( int i=1; i<64; ++i )
    {
      terminateBranch[i] = false;
    }
    groupsPrinted.clear();
  }

  // store the local allocations for each wrapper
  std::vector< size_t > localAllocations;

  // use the first index for the summation of all Wrappers in a Group
  localAllocations.emplace_back( 0 );
  for( auto & view : group.wrappers() )
  {
    size_t const bytesAllocated = view.second->bytesAllocated();
    localAllocations.emplace_back( bytesAllocated );
    localAllocations[0] += bytesAllocated;
  }

  int const numRanks = MpiWrapper::commSize();
  int const numValues = localAllocations.size();

  // storage for the gathered values of localAllocation
  std::vector< size_t > globalAllocations;
  if( MpiWrapper::commRank()==0 )
  {
    globalAllocations.resize( numRanks * numValues );
  }

  MpiWrapper::gather( localAllocations.data(),
                      numValues,
                      globalAllocations.data(),
                      numValues,
                      0,
                      MPI_COMM_GEOS );

  // reduce data across ranks (min, max, sum)
  if( MpiWrapper::commRank()==0 )
  {
    array2d< size_t > allocationReductions( numValues, 3 );
    for( int a=0; a<numValues; ++a )
    {
      allocationReductions( a, 0 ) = std::numeric_limits< size_t >::max();
      allocationReductions( a, 1 ) = 0;
      allocationReductions( a, 2 ) = 0;
      for( int b=0; b<numRanks; ++b )
      {
        int const recvIndex = a + b * numValues;
        size_t const value = globalAllocations[recvIndex];
        allocationReductions( a, 0 ) = std::min( allocationReductions( a, 0 ), value );
        allocationReductions( a, 1 ) = std::max( allocationReductions( a, 1 ), value );
        allocationReductions( a, 2 ) += value;
      }
    }

    // scope for the output of groups
    {
      size_t const * const groupAllocations = allocationReductions[0];

      if( indent==0 )
      {
        //                          1         2         3         4         5         6         7         8         9        10        11
        //        12
        //                 123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        GEOS_LOG_RANK_0( "************************************************************************************************************************" );
        GEOS_LOG_RANK_0( " Data repository memory allocations " );
        GEOS_LOG_RANK_0( "                                                                                       min          max          sum    " );
        GEOS_LOG_RANK_0( "                                                                                   -----------  -----------  -----------" );
      }


      string outputLine;
      int indentChars = 0;
      for( int i=0; i<indent; ++i )
      {
        outputLine += terminateBranch[i]==true ? "   " : "|  ";
        indentChars += 3;
      }
      // put a indent between groups in the tree
      GEOS_LOG_RANK_0( outputLine.c_str()<<"|" );
      indentChars += 1;
      // only allocation data if it is above the threshold...
      if( groupAllocations[0] >= threshold )
      {
        indentChars += 3;
        outputLine += "|--{:.<" + std::to_string( 83-indentChars ) + "} {:>9s}    {:>9s}    {:>9s}";
        GEOS_LOG_RANK_0( GEOS_FMT( outputLine.c_str(),
                                   "[" + group.getName() + "]",
                                   stringutilities::toMetricPrefixString( groupAllocations[0] ) + 'B',
                                   stringutilities::toMetricPrefixString( groupAllocations[1] ) + 'B',
                                   stringutilities::toMetricPrefixString( groupAllocations[2] ) + 'B' ) );
      }
      else  // ...but we still need to output the group name to have a valid tree.
      {
        outputLine += "|--[{:<}]";
        GEOS_LOG_RANK_0( GEOS_FMT( outputLine.c_str(),
                                   group.getName() ) );
      }
    }

    // output the views
    {
      localIndex viewCount = 1;
      for( auto & view : group.wrappers() )
      {
        if( allocationReductions( viewCount, 0 ) > threshold )
        {
          string outputLine;
          int indentChars = 0;
          for( int i=0; i<=indent; ++i )
          {
            outputLine += terminateBranch[i]==true ? "   " : "|  ";
            indentChars += 3;
          }
          indentChars += 5;
          outputLine += "| - {:.<" + std::to_string( 83-indentChars ) + "} {:>9s}    {:>9s}    {:>9s}";
          GEOS_LOG_RANK_0( GEOS_FMT( outputLine.c_str(),
                                     view.second->getName(),
                                     stringutilities::toMetricPrefixString( allocationReductions( viewCount, 0 ) ) + 'B',
                                     stringutilities::toMetricPrefixString( allocationReductions( viewCount, 1 ) ) + 'B',
                                     stringutilities::toMetricPrefixString( allocationReductions( viewCount, 2 ) ) + 'B' ) );
        }
        ++viewCount;
      }
    }
  }

  // process subgroups
  localIndex const numSubGroups = group.numSubGroups();
  localIndex groupCounter = 0;
  for( auto & subGroup : group.getSubGroups() )
  {
    // set flag to indicate whether or not this is the last subgroup so that
    // the ascii art tree may be constructed properly
    terminateBranch[indent+1] = ++groupCounter==numSubGroups ? true : false;

    // check to see that the group hasen't been printed yet.
    if( groupsPrinted.count( subGroup.second ) == 0 )
    {
      // insert the address of the Group into the map that holds printed groups.
      groupsPrinted.insert( subGroup.second );
      printMemoryAllocation( *(subGroup.second), indent + 1, threshold );
    }
  }
  if( indent==0 )
  {
    //         1         2         3         4         5         6         7         8         9        10
    //1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    GEOS_LOG_RANK_0( "****************************************************************************************************" );
  }
}

}
}
