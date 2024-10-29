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

/**
 * @file ConformingVirtualElementOrder1_impl.hpp
 */

#include "mesh/utilities/ComputationalGeometry.hpp"

namespace geos
{
namespace finiteElement
{
/**
 * \private
 */
template< localIndex MCN, localIndex MFN >
template< typename SUBREGION_TYPE >
GEOS_HOST_DEVICE
void ConformingVirtualElementOrder1< MCN, MFN >::
computeProjectors( localIndex const & cellIndex,
                   InputNodeCoords const & nodesCoords,
                   InputCellToNodeMap< SUBREGION_TYPE > const & cellToNodeMap,
                   InputCellToFaceMap< SUBREGION_TYPE > const & elementToFaceMap,
                   InputFaceToNodeMap const & faceToNodeMap,
                   InputFaceToEdgeMap const & faceToEdgeMap,
                   InputEdgeToNodeMap const & edgeToNodeMap,
                   arrayView2d< real64 const > const faceCenters,
                   arrayView2d< real64 const > const faceNormals,
                   arrayView1d< real64 const > const faceAreas,
                   real64 const (&cellCenter)[3],
                   real64 const & cellVolume,
                   localIndex & numSupportPoints,
                   real64 & quadratureWeight,
                   real64 (& basisFunctionsIntegralMean)[MCN],
                   real64 (& stabilizationMatrix)[MCN][MCN],
                   real64 (& basisDerivativesIntegralMean)[MCN][3]
                   )
{
  localIndex const numCellFaces = elementToFaceMap[cellIndex].size();
  localIndex const numCellPoints = cellToNodeMap[cellIndex].size();
  numSupportPoints = numCellPoints;
  quadratureWeight = cellVolume;

  // Compute cell diameter.
  real64 cellDiameter = ConformingVirtualElementOrder1< MCN, MFN >::
                        computeDiameter< 3 >
                          ( nodesCoords, cellToNodeMap[cellIndex], numCellPoints );
  real64 const invCellDiameter = 1.0/cellDiameter;

  // Compute basis functions and scaled monomials integrals on the boundary.
  real64 basisBoundaryIntegrals[ maxSupportPoints ];
  real64 basisTimesNormalBoundaryInt[ maxSupportPoints ][ 3 ];
  real64 basisTimesMonomNormalDerBoundaryInt[ maxSupportPoints ][ 3 ];
  real64 monomBoundaryIntegrals[ 4 ] = { 0.0 };
  // - initialize vectors
  for( localIndex numBasisFunction = 0; numBasisFunction < numCellPoints; ++numBasisFunction )
  {
    basisBoundaryIntegrals[ numBasisFunction ] = 0.0;
    for( localIndex i = 0; i < 3; ++i )
    {
      basisTimesNormalBoundaryInt[ numBasisFunction ][ i ] = 0.0;
    }
    for( localIndex i = 0; i < 3; ++i )
    {
      basisTimesMonomNormalDerBoundaryInt[ numBasisFunction ][ i ] = 0.0;
    }
  }
  // - loop over faces and perform computations on the boundary
  for( localIndex numFace = 0; numFace < numCellFaces; ++numFace )
  {
    localIndex const faceIndex = elementToFaceMap[ cellIndex ][ numFace ];
    real64 const faceArea = faceAreas[ faceIndex ];
    localIndex const numFaceNodes = faceToNodeMap[ faceIndex ].size();
    localIndex faceToNodes[ MFN ];
    localIndex faceToEdges[ MFN ];
    for( localIndex i = 0; i < numFaceNodes; ++i )
    {
      faceToNodes[i] = faceToNodeMap[ faceIndex ][ i ];
      faceToEdges[i] = faceToEdgeMap[ faceIndex ][ i ];
    }
    // - get outward face normal and center
    real64 faceNormal[3] = { faceNormals[faceIndex][0],
                             faceNormals[faceIndex][1],
                             faceNormals[faceIndex][2] };
    real64 const faceCenter[3] { faceCenters[faceIndex][0],
                                 faceCenters[faceIndex][1],
                                 faceCenters[faceIndex][2] };
    // - compute integrals calling auxiliary method
    real64 faceBasisIntegrals[MFN];
    real64 threeDMonomialIntegrals[3] = { 0.0 };
    computeFaceIntegrals( nodesCoords,
                          faceToNodes,
                          faceToEdges,
                          numFaceNodes,
                          faceArea,
                          faceCenter,
                          faceNormal,
                          edgeToNodeMap,
                          invCellDiameter,
                          cellCenter,
                          faceBasisIntegrals,
                          threeDMonomialIntegrals );

    real64 signTestVector[3];
    signTestVector[0] = faceCenters[faceIndex][0] - cellCenter[0];
    signTestVector[1] = faceCenters[faceIndex][1] - cellCenter[1];
    signTestVector[2] = faceCenters[faceIndex][2] - cellCenter[2];
    if( signTestVector[0]*faceNormal[0] +
        signTestVector[1]*faceNormal[1] +
        signTestVector[2]*faceNormal[2] < 0 )
    {
      faceNormal[0] = -faceNormal[0];
      faceNormal[1] = -faceNormal[1];
      faceNormal[2] = -faceNormal[2];
    }
    // - add contributions to integrals of monomials
    monomBoundaryIntegrals[0] += faceArea;
    for( localIndex monomInd = 1; monomInd < 4; ++monomInd )
    {
      monomBoundaryIntegrals[monomInd] += threeDMonomialIntegrals[monomInd-1];
    }
    // - add contributions to integrals of basis functions
    for( localIndex numFaceBasisFunction = 0; numFaceBasisFunction < numFaceNodes;
         ++numFaceBasisFunction )
    {
      localIndex basisFunctionIndex = 0;
      // find the position of the current face vertex within cell vertices
      while( cellToNodeMap[cellIndex][basisFunctionIndex] != faceToNodes[numFaceBasisFunction] )
      {
        ++basisFunctionIndex;
      }

      basisBoundaryIntegrals[basisFunctionIndex] += faceBasisIntegrals[numFaceBasisFunction];
      for( localIndex pos = 0; pos < 3; ++pos )
      {
        basisTimesNormalBoundaryInt[basisFunctionIndex][pos] +=
          faceNormal[pos]*faceBasisIntegrals[numFaceBasisFunction];
        basisTimesMonomNormalDerBoundaryInt[basisFunctionIndex][pos] +=
          faceNormal[pos]*faceBasisIntegrals[numFaceBasisFunction]*invCellDiameter;
      }
    }
  }

  // Compute non constant scaled monomials' integrals on the polyhedron.
  real64 monomInternalIntegrals[3] = { 0.0 };
  for( localIndex numFace = 0; numFace < numCellFaces; ++numFace )
  {
    localIndex const faceIndex = elementToFaceMap[cellIndex][numFace];
    localIndex const numFaceVertices = faceToNodeMap[faceIndex].size();
    localIndex faceToNodes[ MFN ];
    for( localIndex i = 0; i < numFaceVertices; ++i )
    {
      faceToNodes[i] = faceToNodeMap[ faceIndex ][ i ];
    }
    real64 const faceCenter[3] { faceCenters[faceIndex][0],
                                 faceCenters[faceIndex][1],
                                 faceCenters[faceIndex][2] };
    for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
    {
      localIndex numNextVertex = (numVertex+1)%numFaceVertices;
      // compute value of 3D monomials at the quadrature point on the sub-tetrahedron (the
      // barycenter).
      // The result is ((v0 + v1 + faceCenter + cellCenter)/4 - cellCenter) / cellDiameter =
      // = (v0 + v1 + faceCenter - 3*cellcenter)/(4*cellDiameter).
      real64 monomialValues[3];
      for( localIndex pos = 0; pos < 3; ++pos )
      {
        monomialValues[pos] = (nodesCoords( faceToNodes[ numVertex ], pos ) +
                               nodesCoords( faceToNodes[ numNextVertex ], pos ) +
                               faceCenter[ pos ] - 3*cellCenter[ pos ])*invCellDiameter/4.0;
      }
      // compute quadrature weight (the volume of the sub-tetrahedron).
      real64 edgeTangentsMatrix[3][3];
      for( localIndex i = 0; i < 3; ++i )
      {
        edgeTangentsMatrix[0][i] = faceCenter[ i ] - cellCenter[ i ];
        edgeTangentsMatrix[1][i] = nodesCoords( faceToNodes[ numVertex ], i ) - cellCenter[ i ];
        edgeTangentsMatrix[2][i] = nodesCoords( faceToNodes[ numNextVertex ], i ) -
                                   cellCenter[ i ];
      }
      real64 subTetVolume = LvArray::math::abs
                              ( edgeTangentsMatrix[0][0] *
                              ( edgeTangentsMatrix[1][1] * edgeTangentsMatrix[2][2] -
                                edgeTangentsMatrix[1][2] * edgeTangentsMatrix[2][1] ) +
                              edgeTangentsMatrix[1][0] *
                              ( edgeTangentsMatrix[0][2] * edgeTangentsMatrix[2][1] -
                                edgeTangentsMatrix[0][1] * edgeTangentsMatrix[2][2] ) +
                              edgeTangentsMatrix[2][0] *
                              ( edgeTangentsMatrix[0][1] * edgeTangentsMatrix[1][2] -
                                edgeTangentsMatrix[0][2] * edgeTangentsMatrix[1][1] )
                              ) / 6.0;
      for( localIndex i = 0; i < 3; ++i )
      {
        monomInternalIntegrals[i] += monomialValues[i]*subTetVolume;
      }
    }
  }

  // Compute integral mean of basis functions and of derivatives of basis functions.
  // Compute VEM degrees of freedom of the piNabla projection minus the identity (used for
  // stabilizationMatrix).
  real64 const invCellVolume = 1.0/cellVolume;
  real64 const monomialDerivativeInverse = cellDiameter*cellDiameter*invCellVolume;
  real64 piNablaVemDofsMinusIdentity[ MCN ][ MCN ];
  // - compute values of scaled monomials at the vertices (used for piNablaVemDofs)
  real64 monomialVemDofs[ 3 ][ MCN ];
  for( localIndex numVertex = 0; numVertex < numCellPoints; ++numVertex )
  {
    for( localIndex pos = 0; pos < 3; ++pos )
    {
      monomialVemDofs[ pos ][ numVertex ] = invCellDiameter*
                                            (nodesCoords( cellToNodeMap( cellIndex, numVertex ), pos )
                                             - cellCenter[ pos ]);
    }
  }
  for( localIndex numBasisFunction = 0; numBasisFunction < numCellPoints; ++numBasisFunction )
  {
    // - solve linear system to obtain dofs of piNabla proj wrt monomial basis
    real64 piNablaDofs[4];
    piNablaDofs[1] = monomialDerivativeInverse *
                     basisTimesMonomNormalDerBoundaryInt[numBasisFunction][0];
    piNablaDofs[2] = monomialDerivativeInverse *
                     basisTimesMonomNormalDerBoundaryInt[numBasisFunction][1];
    piNablaDofs[3] = monomialDerivativeInverse *
                     basisTimesMonomNormalDerBoundaryInt[numBasisFunction][2];
    piNablaDofs[0] = (basisBoundaryIntegrals[numBasisFunction] -
                      piNablaDofs[1]*monomBoundaryIntegrals[1] -
                      piNablaDofs[2]*monomBoundaryIntegrals[2] -
                      piNablaDofs[3]*monomBoundaryIntegrals[3] )/monomBoundaryIntegrals[0];
    // - integrate piNabla proj and compute integral means
    basisFunctionsIntegralMean[numBasisFunction] =
      piNablaDofs[0] + invCellVolume * (piNablaDofs[1]*monomInternalIntegrals[0]
                                        + piNablaDofs[2]*monomInternalIntegrals[1]
                                        + piNablaDofs[3] * monomInternalIntegrals[2]);
    // - compute integral means of derivatives
    for( localIndex i = 0; i < 3; ++i )
    {
      basisDerivativesIntegralMean[numBasisFunction][i] =
        invCellVolume * basisTimesNormalBoundaryInt[numBasisFunction][i];
    }
    // - compute VEM dofs of piNabla projection
    for( localIndex numVertex = 0; numVertex < numCellPoints; ++numVertex )
    {
      piNablaVemDofsMinusIdentity[ numVertex ][ numBasisFunction ] =
        piNablaDofs[0] + piNablaDofs[1]*monomialVemDofs[ 0 ][ numVertex ] +
        piNablaDofs[2]*monomialVemDofs[ 1 ][ numVertex ] +
        piNablaDofs[3]*monomialVemDofs[ 2 ][ numVertex ];
    }
    piNablaVemDofsMinusIdentity[ numBasisFunction ][ numBasisFunction ] -= 1;
  }

  // Compute stabilization matrix.
  // The result is piNablaVemDofsMinusIdentity^T * piNablaVemDofsMinusIdentity.
  for( localIndex i = 0; i < numCellPoints; ++i )
  {
    for( localIndex j = 0; j < numCellPoints; ++j )
    {
      real64 rowColProd = 0;
      for( localIndex k = 0; k < numCellPoints; ++k )
      {
        rowColProd += piNablaVemDofsMinusIdentity[ k ][ i ]*piNablaVemDofsMinusIdentity[ k ][ j ];
      }
      stabilizationMatrix[ i ][ j ] = cellDiameter*rowColProd;
    }
  }
}

template< localIndex MCN, localIndex MFN >
GEOS_HOST_DEVICE
void ConformingVirtualElementOrder1< MCN, MFN >::
computeFaceIntegrals( InputNodeCoords const & nodesCoords,
                      localIndex const (&faceToNodes)[MFN],
                      localIndex const (&faceToEdges)[MFN],
                      localIndex const & numFaceVertices,
                      real64 const & faceArea,
                      real64 const (&faceCenter)[3],
                      real64 const (&faceNormal)[3],
                      InputEdgeToNodeMap const & edgeToNodes,
                      real64 const & invCellDiameter,
                      real64 const (&cellCenter)[3],
                      real64 (& basisIntegrals)[MFN],
                      real64 (& threeDMonomialIntegrals)[3] )
{
  // Rotate the face.
  //  - compute rotation matrix.
  real64 faceRotationMatrix[ 3 ][ 3 ];
  computationalGeometry::RotationMatrix_3D( faceNormal, faceRotationMatrix );
  //  - below we compute the diameter, the rotated vertices and the rotated center.
  real64 faceRotatedVertices[ MFN ][ 2 ];
  real64 faceDiameter = 0;
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    // apply the transpose (that is the inverse) of the rotation matrix to face vertices.
    // NOTE:
    // the second and third rows of the transpose of the rotation matrix rotate on the 2D face.
    faceRotatedVertices[numVertex][0] =
      faceRotationMatrix[ 0 ][ 1 ]*nodesCoords( faceToNodes[ numVertex ], 0 ) +
      faceRotationMatrix[ 1 ][ 1 ]*nodesCoords( faceToNodes[ numVertex ], 1 ) +
      faceRotationMatrix[ 2 ][ 1 ]*nodesCoords( faceToNodes[ numVertex ], 2 );
    faceRotatedVertices[numVertex][1] =
      faceRotationMatrix[ 0 ][ 2 ]*nodesCoords( faceToNodes[ numVertex ], 0 ) +
      faceRotationMatrix[ 1 ][ 2 ]*nodesCoords( faceToNodes[ numVertex ], 1 ) +
      faceRotationMatrix[ 2 ][ 2 ]*nodesCoords( faceToNodes[ numVertex ], 2 );
  }
  faceDiameter = ConformingVirtualElementOrder1< MCN, MFN >::
                 computeDiameter< 2 >( faceRotatedVertices,
                                       numFaceVertices );
  real64 const invFaceDiameter = 1.0/faceDiameter;
  // - rotate the face centroid as done for the vertices.
  real64 faceRotatedCentroid[2];
  faceRotatedCentroid[0] =
    faceRotationMatrix[ 0 ][ 1 ]*faceCenter[0] +
    faceRotationMatrix[ 1 ][ 1 ]*faceCenter[1] +
    faceRotationMatrix[ 2 ][ 1 ]*faceCenter[2];
  faceRotatedCentroid[1] =
    faceRotationMatrix[ 0 ][ 2 ]*faceCenter[0] +
    faceRotationMatrix[ 1 ][ 2 ]*faceCenter[1] +
    faceRotationMatrix[ 2 ][ 2 ]*faceCenter[2];
  // - compute edges' lengths, outward pointing normals and local edge-to-nodes map.
  real64 edgeOutwardNormals[ MFN ][ 2 ];
  real64 edgeLengths[ MFN ];
  localIndex localEdgeToNodes[ MFN ][ 2 ];
  for( localIndex numEdge = 0; numEdge < numFaceVertices; ++numEdge )
  {
    if( edgeToNodes( faceToEdges[numEdge], 0 ) == faceToNodes[ numEdge ] )
    {
      localEdgeToNodes[ numEdge ][ 0 ] = numEdge;
      localEdgeToNodes[ numEdge ][ 1 ] = (numEdge+1)%numFaceVertices;
    }
    else
    {
      localEdgeToNodes[ numEdge ][ 0 ] = (numEdge+1)%numFaceVertices;
      localEdgeToNodes[ numEdge ][ 1 ] = numEdge;
    }
    real64 edgeTangent[2];
    edgeTangent[0] = faceRotatedVertices[(numEdge+1)%numFaceVertices][0] -
                     faceRotatedVertices[numEdge][0];
    edgeTangent[1] = faceRotatedVertices[(numEdge+1)%numFaceVertices][1] -
                     faceRotatedVertices[numEdge][1];
    edgeOutwardNormals[numEdge][0] = edgeTangent[1];
    edgeOutwardNormals[numEdge][1] = -edgeTangent[0];
    real64 signTestVector[2];
    signTestVector[0] = faceRotatedVertices[numEdge][0] - faceRotatedCentroid[0];
    signTestVector[1] = faceRotatedVertices[numEdge][1] - faceRotatedCentroid[1];
    if( signTestVector[0]*edgeOutwardNormals[numEdge][0] +
        signTestVector[1]*edgeOutwardNormals[numEdge][1] < 0 )
    {
      edgeOutwardNormals[numEdge][0] = -edgeOutwardNormals[numEdge][0];
      edgeOutwardNormals[numEdge][1] = -edgeOutwardNormals[numEdge][1];
    }
    edgeLengths[numEdge] = LvArray::math::sqrt< real64 >( edgeTangent[0]*edgeTangent[0] +
                                                          edgeTangent[1]*edgeTangent[1] );
    edgeOutwardNormals[numEdge][0] /= edgeLengths[numEdge];
    edgeOutwardNormals[numEdge][1] /= edgeLengths[numEdge];
  }

  // Compute boundary quadrature weights (also equal to the integrals of basis functions on the
  // boundary).
  real64 boundaryQuadratureWeights[ MFN ];
  for( localIndex numWeight = 0; numWeight < numFaceVertices; ++numWeight )
    boundaryQuadratureWeights[numWeight] = 0.0;
  for( localIndex numEdge = 0; numEdge < numFaceVertices; ++numEdge )
  {
    boundaryQuadratureWeights[ localEdgeToNodes[ numEdge ][ 0 ] ] += 0.5*edgeLengths[numEdge];
    boundaryQuadratureWeights[ localEdgeToNodes[ numEdge ][ 1 ] ] += 0.5*edgeLengths[numEdge];
  }

  // Compute scaled monomials' integrals on edges.
  real64 monomBoundaryIntegrals[3] = { 0.0 };
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    monomBoundaryIntegrals[0] += boundaryQuadratureWeights[ numVertex ];
    monomBoundaryIntegrals[1] += (faceRotatedVertices[ numVertex ][ 0 ] - faceRotatedCentroid[0]) *
                                 invFaceDiameter*boundaryQuadratureWeights[ numVertex ];
    monomBoundaryIntegrals[2] += (faceRotatedVertices[ numVertex ][ 1 ] - faceRotatedCentroid[1]) *
                                 invFaceDiameter*boundaryQuadratureWeights[ numVertex ];
  }

  // Compute non constant 2D and 3D scaled monomials' integrals on the face.
  real64 monomInternalIntegrals[2] = { 0.0 };
  for( localIndex numSubTriangle = 0; numSubTriangle < numFaceVertices; ++numSubTriangle )
  {
    localIndex const nextVertex = (numSubTriangle+1)%numFaceVertices;
    // - compute value of 2D monomials at the quadrature point on the sub-triangle (the
    //   barycenter).
    //   The result is ((v(0)+v(1)+faceCenter)/3 - faceCenter) / faceDiameter =
    //   = (v(0) + v(1) - 2*faceCenter)/(3*faceDiameter).
    real64 monomialValues[2];
    for( localIndex i = 0; i < 2; ++i )
    {
      monomialValues[i] = (faceRotatedVertices[numSubTriangle][i] +
                           faceRotatedVertices[nextVertex][i] -
                           2.0*faceRotatedCentroid[i]) / (3.0*faceDiameter);
    }
    // compute value of 3D monomials at the quadrature point on the sub-triangle (the
    // barycenter).  The result is
    // ((v(0) + v(1) + faceCenter)/3 - cellCenter)/cellDiameter.
    real64 threeDMonomialValues[3];
    for( localIndex i = 0; i < 3; ++i )
    {
      threeDMonomialValues[i] = ( (faceCenter[i] +
                                   nodesCoords[faceToNodes[ numSubTriangle ]][i] +
                                   nodesCoords[faceToNodes[ nextVertex ]][i]) / 3.0 -
                                  cellCenter[i] ) * invCellDiameter;
    }
    // compute quadrature weight associated to the quadrature point (the area of the
    // sub-triangle).
    real64 edgesTangents[2][2];               // used to compute the area of the sub-triangle
    for( localIndex i = 0; i < 2; ++i )
    {
      edgesTangents[0][i] = faceRotatedVertices[numSubTriangle][i] - faceRotatedCentroid[i];
    }
    for( localIndex i = 0; i < 2; ++i )
    {
      edgesTangents[1][i] = faceRotatedVertices[nextVertex][i] - faceRotatedCentroid[i];
    }
    real64 subTriangleArea = 0.5*LvArray::math::abs
                               ( edgesTangents[0][0]*edgesTangents[1][1] -
                               edgesTangents[0][1]*edgesTangents[1][0] );
    // compute the integrals on the sub-triangle and add it to the global integrals
    for( localIndex i = 0; i < 2; ++i )
    {
      monomInternalIntegrals[ i ] += monomialValues[ i ]*subTriangleArea;
    }
    for( localIndex i = 0; i < 3; ++i )
    {
      // threeDMonomialIntegrals is assumed to be initialized to 0 by the caller
      threeDMonomialIntegrals[ i ] += threeDMonomialValues[ i ]*subTriangleArea;
    }
  }

  // Compute integral of basis functions times normal derivative of monomials on the boundary.
  real64 basisTimesMonomNormalDerBoundaryInt[ MFN ][ 2 ];
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    for( localIndex i = 0; i < 2; ++i )
    {
      basisTimesMonomNormalDerBoundaryInt[ numVertex ][ i ] = 0.0;
    }
  }
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    for( localIndex i = 0; i < 2; ++i )
    {
      real64 thisEdgeIntTimesNormal_i = edgeOutwardNormals[numVertex][i]*edgeLengths[numVertex];
      basisTimesMonomNormalDerBoundaryInt[ localEdgeToNodes[ numVertex ][ 0 ] ][i] += thisEdgeIntTimesNormal_i;
      basisTimesMonomNormalDerBoundaryInt[ localEdgeToNodes[ numVertex ][ 1 ] ][i] += thisEdgeIntTimesNormal_i;
    }
  }
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    for( localIndex i = 0; i < 2; ++i )
    {
      basisTimesMonomNormalDerBoundaryInt[ numVertex ][ i ] *= 0.5*invFaceDiameter;
    }
  }

  // Compute integral mean of basis functions on this face.
  real64 const invFaceArea = 1.0/faceArea;
  real64 const monomialDerivativeInverse = (faceDiameter*faceDiameter)*invFaceArea;
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    real64 piNablaDofs[ 3 ];
    piNablaDofs[ 1 ] = monomialDerivativeInverse *
                       basisTimesMonomNormalDerBoundaryInt[ numVertex ][ 0 ];
    piNablaDofs[ 2 ] = monomialDerivativeInverse *
                       basisTimesMonomNormalDerBoundaryInt[ numVertex ][ 1 ];
    piNablaDofs[ 0 ] = (boundaryQuadratureWeights[ numVertex ] -
                        piNablaDofs[ 1 ]*monomBoundaryIntegrals[ 1 ] -
                        piNablaDofs[ 2 ]*monomBoundaryIntegrals[ 2 ])/monomBoundaryIntegrals[ 0 ];
    basisIntegrals[ numVertex ] = piNablaDofs[ 0 ]*faceArea +
                                  (piNablaDofs[ 1 ]*monomInternalIntegrals[ 0 ] +
                                   piNablaDofs[ 2 ]*monomInternalIntegrals[ 1 ]);
  }
}
}
}
