/*
 * BoundaryConditionManager.hpp
 *
 *  Created on: May 26, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_
#include "common/DataTypes.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "BoundaryConditionBase.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const boundaryConditionMananger("BoundaryConditionMananger");
}
}

class BoundaryConditionManager : public dataRepository::ManagedGroup
{
public:
  BoundaryConditionManager( string const & name, dataRepository::ManagedGroup * const parent );
  virtual ~BoundaryConditionManager();

  static BoundaryConditionManager & get();

  void ReadXMLsub( xmlWrapper::xmlNode const & targetNode );


  template< typename BCFunctionPtr, typename ... ARGS>
  void ApplyBoundaryCondition( BCFunctionPtr boundaryConditionFunctionPtr,
                               dataRepository::ManagedGroup & object,
                               std::string const & fieldName,
                               real64 const time,
                               ARGS & ... args );

  template< typename Solver, typename BCFunctionPtr, typename ... ARGS>
  void ApplyBoundaryCondition( Solver* solverPtr,
                               BCFunctionPtr boundaryConditionFunctionPtr,
                               dataRepository::ManagedGroup & object,
                               std::string const & fieldName,
                               real64 const time,
                               ARGS ... args );
};


template< typename BCFunctionPtr, typename ... ARGS>
void BoundaryConditionManager::ApplyBoundaryCondition( BCFunctionPtr boundaryConditionFunctionPtr,
                                                       dataRepository::ManagedGroup & object,
                                                       std::string const & fieldName,
                                                       real64 const time,
                                                       ARGS & ... args )
{
  dataRepository::ManagedGroup const & sets = object.GetGroup(dataRepository::keys::sets);

  // iterate over all boundary conditions.
  forSubGroups<BoundaryConditionBase>( [&]( BoundaryConditionBase & bc ) -> void
  {
    if( time >= bc.GetStartTime() && time < bc.GetEndTime() && ( bc.GetFieldName()==fieldName) )
    {
      string_array setNames = bc.GetSetNames();
      for( auto & setName : setNames )
      {
        dataRepository::ViewWrapper<lSet> * setWrapper = getWrapperPtr<lSet>(setName);
        if( setWrapper != nullptr )
        {
          lSet const & set = setWrapper->reference();
          (*boundaryConditionFunctionPtr)( &bc, set, time, args... );
        }
      }
    }
  });
}


template< typename Solver, typename BCFunctionPtr, typename ... ARGS>
void BoundaryConditionManager::ApplyBoundaryCondition( Solver* solverPtr,
                                                       BCFunctionPtr boundaryConditionFunctionPtr,
                                                       dataRepository::ManagedGroup & object,
                                                       std::string const & fieldName,
                                                       real64 const time,
                                                       ARGS... args )
{
  dataRepository::ManagedGroup const & sets = object.GetGroup("Sets");

  // iterate over all boundary conditions.
  forSubGroups( [&]( BoundaryConditionBase & bc ) -> void
  {
    if( time >= bc.GetStartTime() && time < bc.GetEndTime() && ( bc.GetFieldName()==fieldName) )
    {
      for(localIndex i =0; i < bc.m_setNames.size(); ++i)
      {
        int findSet = 1;

        string_array setNames = bc.GetSetNames();
        for( auto & setName : setNames )
        {
          lSet const * const set = object.hasGroup(setName);
          if( set != nullptr )
          {
            (solverPtr->*boundaryConditionFunctionPtr)(object,bc,*set,time,args...);
          }
        }
      }
    }
  });
}



} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_ */
