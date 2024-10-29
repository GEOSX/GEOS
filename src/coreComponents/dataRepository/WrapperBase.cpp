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

/** @file */

#include "WrapperBase.hpp"

#include "Group.hpp"
#include "RestartFlags.hpp"
#include "WrapperContext.hpp"


namespace geos
{
namespace dataRepository
{


WrapperBase::WrapperBase( string const & name,
                          Group & parent,
                          string const & rtTypeName ):
  m_name( name ),
  m_parent( &parent ),
  m_sizedFromParent( 1 ),
  m_restart_flags( RestartFlags::WRITE_AND_READ ),
  m_plotLevel( PlotLevel::NOPLOT ),
  m_inputFlag( InputFlags::INVALID ),
  m_successfulReadFromInput( false ),
  m_description(),
  m_rtTypeName( rtTypeName ),
  m_registeringObjects(),
  m_conduitNode( parent.getConduitNode()[ name ] ),
  m_dataContext( std::make_unique< WrapperContext >( *this ) )
{}


WrapperBase::~WrapperBase()
{}

void WrapperBase::resize()
{
  resize( m_parent->size());
}

void WrapperBase::copyWrapperAttributes( WrapperBase const & source )
{
  m_sizedFromParent = source.m_sizedFromParent;
  m_restart_flags = source.m_restart_flags;
  m_plotLevel  = source.m_plotLevel;
  m_inputFlag = source.m_inputFlag;
  m_description = source.m_description;
  m_rtTypeName = source.m_rtTypeName;
}

string WrapperBase::getPath() const
{
  // In the Conduit node hierarchy everything begins with 'Problem', we should change it so that
  // the ProblemManager actually uses the root Conduit Node but that will require a full rebaseline.
  string const noProblem = m_conduitNode.path().substr( std::strlen( dataRepository::keys::ProblemManager ) - 1 );
  return noProblem.empty() ? "/" : noProblem;
}

string WrapperBase::dumpInputOptions( bool const outputHeader ) const
{
  string rval;
  if( outputHeader )
  {
    rval.append( "  |         name         |  opt/req  | Description \n" );
    rval.append( "  |----------------------|-----------|-----------------------------------------\n" );
  }

  if( getInputFlag() == InputFlags::OPTIONAL || getInputFlag() == InputFlags::REQUIRED )
  {
    rval.append( GEOS_FMT( "  | {:20} | {:9} | {} \n", getName(), InputFlagToString( getInputFlag() ), getDescription() ) );
  }

  return rval;
}

#if defined(USE_TOTALVIEW_OUTPUT)
int WrapperBase::setTotalviewDisplay() const
{
  //std::cout<<"exectuing WrapperBase::setTotalviewDisplay()"<<std::endl;
//  TV_ttf_add_row("TYPE", TV_ttf_type_ascii_string, type.c_str() );
  TV_ttf_add_row( "m_name", LvArray::system::demangle< string >().c_str(), &m_name );
  TV_ttf_add_row( "m_parent", LvArray::system::demangle< Group >().c_str(), m_parent );
  TV_ttf_add_row( "m_sizedFromParent", "int", &m_sizedFromParent );
  TV_ttf_add_row( "m_restart_flags", LvArray::system::demangle< RestartFlags >().c_str(), &m_restart_flags );
  TV_ttf_add_row( "m_plotLevel", LvArray::system::demangle< PlotLevel >().c_str(), &m_plotLevel );
  TV_ttf_add_row( "m_inputFlag", LvArray::system::demangle< InputFlags >().c_str(), &m_inputFlag );
  TV_ttf_add_row( "m_description", LvArray::system::demangle< string >().c_str(), &m_description );
  size_t junk = m_registeringObjects.size();
  TV_ttf_add_row( "m_registeringObjects",
                  totalview::format< string, size_t >( 1, &junk ).c_str(),
                  m_registeringObjects.data() );

  return 0;
}
#endif

void WrapperBase::createDataContext( xmlWrapper::xmlNode const & targetNode,
                                     xmlWrapper::xmlNodePos const & nodePos )
{
  xmlWrapper::xmlAttribute att = targetNode.attribute( m_name.c_str() );
  xmlWrapper::xmlAttributePos attPos = nodePos.getAttributeLine( m_name );
  if( nodePos.isFound() && attPos.isFound() && !att.empty() )
  {
    m_dataContext = std::make_unique< DataFileContext >( targetNode, att, attPos );
  }
}


}
} /* namespace geos */

#if defined(USE_TOTALVIEW_OUTPUT)
/**
 * @brief Global function correlated with WrapperBase to be called by Totalview when displaying
 *        a WrapperBase as a VieWrapper<T>
 * @param wrapper A pointer to the wrapper that will be displayed.
 * @return 0
 */
int TV_ttf_display_type( const geos::dataRepository::WrapperBase * wrapper )
{
  if( wrapper!=nullptr )
  {
    //std::cout<<"displaying WrapperBase "<<wrapper->getName()<<" as "<<wrapper->totalviewTypeName()<<std::endl;
// keep this and try to make it work later on.
//    rval = TV_ttf_add_row( "casted_this", wrapper->totalviewTypeName().c_str(), wrapper );
    wrapper->setTotalviewDisplay();
  }
  return 0;
}
#endif
