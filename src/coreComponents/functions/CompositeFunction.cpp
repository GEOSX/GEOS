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

#include "FunctionManager.hpp"
#include "CompositeFunction.hpp"

namespace geos
{

namespace dataRepository
{
namespace keys
{
string const functionNames = "functionNames";
string const variableNames = "variableNames";
string const expression = "expression";
}
}

using namespace dataRepository;


CompositeFunction::CompositeFunction( const string & name,
                                      Group * const parent ):
  FunctionBase( name, parent ),
  parserContext(),
  parserExpression(),
  m_numSubFunctions(),
  m_subFunctions()
{
  registerWrapper( keys::functionNames, &m_functionNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of source functions. The order must match the variableNames argument." );

  registerWrapper( keys::variableNames, &m_variableNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of variables in expression" );

  registerWrapper( keys::expression, &m_expression ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Composite math expression" );
}

CompositeFunction::~CompositeFunction()
{}

void CompositeFunction::initializeFunction()
{
  // Register variables
  for( localIndex ii=0; ii<m_variableNames.size(); ++ii )
  {
    parserContext.addVariable( m_variableNames[ii].c_str(), static_cast< int >(ii * sizeof(double)));
  }

  // Add built in constants/functions (PI, E, sin, cos, ceil, exp, etc.),
  // compile
  parserContext.addBuiltIns();
  mathpresso::Error err = parserExpression.compile( parserContext, m_expression.c_str(), mathpresso::kNoOptions );
  GEOS_ERROR_IF( err != mathpresso::kErrorOk, "JIT Compiler Error" );

  // Grab pointers to sub functions
  FunctionManager & functionManager = FunctionManager::getInstance();
  m_numSubFunctions = LvArray::integerConversion< localIndex >( m_functionNames.size());
  for( localIndex ii=0; ii<m_numSubFunctions; ++ii )
  {
    m_subFunctions.emplace_back( &functionManager.getGroup< FunctionBase >( m_functionNames[ii] ) );
  }
}

void CompositeFunction::evaluate( dataRepository::Group const & group,
                                  real64 const time,
                                  SortedArrayView< localIndex const > const & set,
                                  arrayView1d< real64 > const & result ) const
{
  // Evaluate each of the subFunctions independently and place the results into
  // a temporary field
  array1d< real64_array > subFunctionResults;
  for( localIndex ii=0; ii<m_numSubFunctions; ++ii )
  {
    real64_array tmp( result.size());
    m_subFunctions[ii]->evaluate( group, time, set, tmp );
    subFunctionResults.emplace_back( std::move( tmp ));
  }

  // Evaluate the symbolic math
  forAll< serialPolicy >( set.size(), [&, result, set]( localIndex const i )
  {
    localIndex const ii = set[ i ];
    real64 functionResults[m_maxNumSubFunctions];
    for( localIndex jj=0; jj<m_numSubFunctions; ++jj )
    {
      functionResults[jj] = subFunctionResults[jj][ii];
    }
    result[ii] = parserExpression.evaluate( reinterpret_cast< void * >( functionResults ));
  } );
}

real64 CompositeFunction::evaluate( real64 const * const input ) const
{
  real64 functionResults[m_maxNumSubFunctions];

  for( localIndex ii=0; ii<m_numSubFunctions; ++ii )
  {
    functionResults[ii] = m_subFunctions[ii]->evaluate( input );
  }

  return parserExpression.evaluate( reinterpret_cast< void * >( functionResults ));
}

REGISTER_CATALOG_ENTRY( FunctionBase, CompositeFunction, string const &, Group * const )

} // namespace geos
