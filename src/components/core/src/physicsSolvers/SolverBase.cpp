/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "SolverBase.hpp"


namespace geosx
{

using namespace dataRepository;

SolverBase::SolverBase( std::string const & name,
                        ManagedGroup * const parent ):
  ManagedGroup( name, parent ),
  m_verboseLevel(0),
  m_gravityVector( R1Tensor(0.0) )
{
  this->RegisterViewWrapper( viewKeyStruct::verboseLevelString, &m_verboseLevel, 0 );
  this->RegisterViewWrapper( viewKeyStruct::gravityVectorString, &m_gravityVector, 0 );

  if( this->globalGravityVector() != nullptr )
  {
    m_gravityVector=*globalGravityVector();
  }
}

SolverBase::~SolverBase()
{}

SolverBase::CatalogInterface::CatalogType& SolverBase::GetCatalog()
{
  static SolverBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void SolverBase::FillDocumentationNode()
{


  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName(this->CatalogName());    // If this method lived in Managed
                                            // groups, this could be done
                                            // automatically
  docNode->setSchemaType("Node");

  docNode->AllocateChildNode( keys::courant,
                              keys::courant,
                              -1,
                              "real64",
                              "real64",
                              "courant Number",
                              "courant Number",
                              "0.7",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::maxDt,
                              keys::maxDt,
                              -1,
                              "real64",
                              "real64",
                              "Maximum Stable Timestep",
                              "Maximum Stable Timestep",
                              "0.0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::verboseLevelString,
                              viewKeyStruct::verboseLevelString,
                              -1,
                              "integer",
                              "integer",
                              "verbosity level",
                              "verbosity level",
                              "0",
                              "",
                              0,
                              1,
                              0 );


}



void SolverBase::TimeStep( real64 const& time_n,
                           real64 const& dt,
                           const int cycleNumber,
                           ManagedGroup * domain )
{
}



void SolverBase::CreateChild( string const & childKey, string const & childName )
{
  if( CatalogInterface::hasKeyName(childKey) )
  {
    std::cout << "Adding Solver of type " << childKey << ", named " << childName << std::endl;
    this->RegisterGroup( childName, CatalogInterface::Factory( childKey, childName, this ) );
  }
}


R1Tensor const * SolverBase::globalGravityVector() const
{
  R1Tensor const * rval = nullptr;
  if( getParent()->getName() == "Solvers" )
  {
    rval = &(getParent()->
             group_cast<SolverBase const *>()->
             getGravityVector());
  }

  return rval;
}


} /* namespace ANST */
