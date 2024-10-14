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
 * @file LaplaceFEM.cpp
 */

#include "mesh/DomainPartition.hpp"

// Source includes
#include "LaplaceFEM.hpp"
#include "LaplaceFEMKernels.hpp"

namespace geos
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;

/*----------------------------------------------------------------------------------
 * LaplaceFEM: Solving Laplace's partial differential equation with finite elements
 * ---------------------------------------------------------------------------------
 *
 * What does this solver do?
 * --------------------------
 *
 * This solver finds a solution f(x,y,z) to the Laplace equation: div ( grad ( f )) = 0
 * This common elliptic PDE represents the solution of a steady-state heat transfer, for instance.
 *
 * Where can I find an example of what it does?
 * --------------------------------------------
 *
 * Integrated tests associated to this solver are found in the ./integratedTests/ folder
 * These tests consist of computing the steady-state temperature profile in a simple cube-shaped domain
 * with fixed temperatures applied on two opposite cube faces ("Dirichlet" boundary conditions: imposing a value).
 * Feel free to run these tests cases, check out the XML input files, and inspect the output.
 *
 * Implementation: before we start:
 * ---------------------------------
 * In this implementation, the solution function (called above f) is called m_fieldName.
 * The variable m_fieldName is a string that points to a data container (an array) that
 * holds the numerical values of the PDE solution for each location at which f is evaluated.
 *
 * Let's take a look at the implementation step by step.
 *
 * ---------------------------------------------------------------------------------
 */


/* CONSTRUCTOR
   First, let us inspect the constructor of a "LaplaceFEM" object.
   This constructor does three important things:
   1 - It constructs an instance of the LaplaceFEM class (here: using the SolverBase constructor and passing through the arguments).
   2 - It sets some default values for the LaplaceFEM-specific private variables (here: m_fieldName and m_timeIntegrationOption).
   3 - It creates and activates a "registerWrapper" for each private variable.
   This is where the private variables are declared either as REQUIRED or OPTIONAL.
   An error is thrown if a REQUIRED variable is not specified in the XML file,
   along with the description of this variable and possible enum values if relevant.
   The description that is set is used in auto-generated documentation and console error messages.
 */

//START_SPHINX_INCLUDE_CONSTRUCTOR
LaplaceFEM::LaplaceFEM( const string & name,
                        Group * const parent ):
  LaplaceBaseH1( name, parent )
{}
//END_SPHINX_INCLUDE_CONSTRUCTOR

LaplaceFEM::~LaplaceFEM()
{
  // TODO Auto-generated destructor stub
}

/* SETUP SYSTEM
   Setting up the system using the base class method
 */
void LaplaceFEM::setupSystem( DomainPartition & domain,
                              DofManager & dofManager,
                              CRSMatrix< real64, globalIndex > & localMatrix,
                              ParallelVector & rhs,
                              ParallelVector & solution,
                              bool const setSparsity )
{
  GEOS_MARK_FUNCTION;
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    string const dofKey = dofManager.getKey( m_fieldName );
    arrayView1d< globalIndex const > const &
    dofIndex = nodeManager.getReference< globalIndex_array >( dofKey );

    SparsityPattern< globalIndex > sparsityPattern( dofManager.numLocalDofs(),
                                                    dofManager.numGlobalDofs(),
                                                    8*8*3 );

    finiteElement::fillSparsity< CellElementSubRegion,
                                 LaplaceFEMKernel >( mesh,
                                                     regionNames,
                                                     this->getDiscretizationName(),
                                                     dofIndex,
                                                     dofManager.rankOffset(),
                                                     sparsityPattern );

    sparsityPattern.compress();
    localMatrix.assimilate< parallelDevicePolicy<> >( std::move( sparsityPattern ) );
  } );
}


/*
   ASSEMBLE SYSTEM
   This is the most important method to assemble the matrices needed before sending them to our solver.
   For a system A.x = B (with x the unknown), here, we use:
   - A : "localMatrix" this represents a Compressed Row Storage (optimized for sparse) matrix of real64 values associated with their index,
   - B : "localRhs" this represents a vector (1d array) of real64 numbers specified at the equation's right-hand side.
   The "local" prefix indicates that we are working on a local problem here, and the parallelization is performed at a higher level.
   This assembly step collects all the information needed to create the matrices localMatrix and localRhs, and the computation of values
   is done in a specific Laplace kernel optimized for parallel performance. Here we:
   1 - identify and point to the mesh of this domain,
   2 - find the node manager of this mesh,
   3 - extract the indices of the nodes that will be solved for (ie. the degrees of freedom or "dof")
   4 - pass all this information to a Laplace-specific finite element computation kernel.
   The call to the kernel is a templated call designed for performance (we will not explain the kernel here).
   See the implementation in LaplaceFEMKernel.cpp.
 */
//START_SPHINX_INCLUDE_ASSEMBLY
void LaplaceFEM::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                 real64 const dt,
                                 DomainPartition & domain,
                                 DofManager const & dofManager,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    string const dofKey = dofManager.getKey( m_fieldName );
    arrayView1d< globalIndex const > const &
    dofIndex =  nodeManager.getReference< array1d< globalIndex > >( dofKey );

    LaplaceFEMKernelFactory kernelFactory( dofIndex, dofManager.rankOffset(), localMatrix, localRhs, dt, m_fieldName );

    string const dummyString = "dummy";
    finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< >,
                                    constitutive::NullModel,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            this->getDiscretizationName(),
                                                            dummyString,
                                                            kernelFactory );

  } );

}
//END_SPHINX_INCLUDE_ASSEMBLY

//START_SPHINX_INCLUDE_REGISTER
REGISTER_CATALOG_ENTRY( SolverBase, LaplaceFEM, string const &, Group * const )
//END_SPHINX_INCLUDE_REGISTER
} /* namespace geos */
