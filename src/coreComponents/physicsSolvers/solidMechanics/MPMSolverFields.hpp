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
 * @file MPMSolverBaseFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERBASEFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERBASEFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace mpm
{

DECLARE_FIELD( particleDeleteFlag,
               "particleDeleteFlag",
               array1d< int >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "An array that remembers particles that should be deleted at the end of the time step." );

DECLARE_FIELD( particleCrystalHealFlag,
               "particleCrystalHealFlag",
               array1d< int >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that remembers particles that are undergoing crystal healing for the mpm event." );

DECLARE_FIELD( particleMaterialType,
               "particleMaterialType",
               array1d< int >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores index of particle material type (assumes single material for each particle region)." );

DECLARE_FIELD( particleMass,
               "particleMass",
               array1d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle masses." );

DECLARE_FIELD( particleWavespeed,
               "particleWavespeed",
               array1d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle wavespeeds." );

DECLARE_FIELD( particleHeatCapacity,
               "particleHeatCapacity",
               array1d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle temperature." );

DECLARE_FIELD( particleReferenceTemperature,
               "particleReferenceTemperature",
               array1d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle reference temperature." );

DECLARE_FIELD( particleInternalEnergy,
               "particleInternalEnergy",
               array1d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle internal energy." );

DECLARE_FIELD( particleKineticEnergy,
               "particleKineticEnergy",
               array1d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle kinetic energy." );

DECLARE_FIELD( particleArtificialViscosity,
               "particleArtificialViscosity",
               array1d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle internal energy." );

DECLARE_FIELD( particleSPHJacobian,
               "particleSPHJacobian",
               array1d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that stores particle SPH computed jacobian." );

DECLARE_FIELD( particleReferenceVolume,
               "particleReferenceVolume",
               array1d< real64 >,
               0.0,
               NOPLOT,
               WRITE_AND_READ,
               "ParticleReferenceVolume" );

DECLARE_FIELD( particleReferencePorosity,
               "particleReferencePorosity",
               array1d< real64 >,
               0.0,
               NOPLOT,
               WRITE_AND_READ,
               "ParticleReferencePorosity" );

DECLARE_FIELD( particleReferenceRVectors,
               "particleReferenceRVectors",
               array3d< real64 >,
               0.0,
               NOPLOT,
               WRITE_AND_READ,
               "ParticleReferenceRVectors" );

DECLARE_FIELD( particleDeformationGradient,
               "particleDeformationGradient",
               array3d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleDeformationGradient" );

DECLARE_FIELD( particleFDot,
               "particleFDot",
               array3d< real64 >,
               0.0,
               NOPLOT,
               WRITE_AND_READ,
               "Material time derivative of the particle deformation gradient." );

DECLARE_FIELD( particleVelocityGradient,
               "particleVelocityGradient",
               array3d< real64 >,
               0.0,
               NOPLOT,
               WRITE_AND_READ,
               "ParticleVelocityGradient" );

DECLARE_FIELD( particleStress,
               "particleStress",
               array2d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle stresses in Voigt notation." );

DECLARE_FIELD( particleBodyForce,
               "particleBodyForce",
               array2d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle body forces." );

DECLARE_FIELD( particlePlasticStrain,
               "particlePlasticStrain",
               array2d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle plastic strain in Voigt notation." );

DECLARE_FIELD( particleDensity,
               "particleDensity",
               array1d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle densities." );

DECLARE_FIELD( particleDamageGradient,
               "particleDamageGradient",
               array2d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle damage gradients as calculated with an SPH kernel." );

DECLARE_FIELD( particleSurfaceFlag,
               "particleSurfaceFlag",
               array1d< int >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "An array that holds particle surface flags." );

DECLARE_FIELD( particleSphF,
               "particleSphF",
               array3d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleSphF" );

DECLARE_FIELD( particleOverlap,
               "particleOverlap",
               array1d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleOverlap" );

DECLARE_FIELD( particleReferencePosition,
               "particleReferencePosition",
               array2d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleReferencePosition" );

DECLARE_FIELD( particleReferenceMaterialDirection, 
               "particleReferenceMaterialDirection", 
               array2d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleReferenceMaterialDirection" );

DECLARE_FIELD( particleReferenceSurfaceNormal, 
               "particleReferenceSurfaceNormal", 
               array2d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleReferenceSurfaceNormal" );

DECLARE_FIELD( particleReferenceSurfacePosition, 
               "particleReferenceSurfacePosition", 
               array2d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleReferenceSurfacePosition" );

DECLARE_FIELD( particleReferenceSurfaceTraction, 
               "particleReferenceSurfaceTraction", 
               array2d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleReferenceSurfaceTraction" );

DECLARE_FIELD( particleCohesiveForce,
               "particleCohesiveForce",
               array2d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleCohesiveForce" );

DECLARE_FIELD( particleCohesiveZoneFlag, 
               "particleCohesiveZoneFlag", 
               array1d< int >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleCohesiveZoneFlag" );

DECLARE_FIELD( particleReferenceMappedNodes, 
               "particleReferenceMappedNodes", 
               array2d< globalIndex >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleReferenceMappedNodes" );

DECLARE_FIELD( particleReferenceShapeFunctionValues, 
               "particleReferenceShapeFunctionValues", 
               array2d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleReferenceShapeFunctionValues" );

DECLARE_FIELD( particleReferenceShapeFunctionGradientValues, 
               "particleReferenceShapeFunctionGradientValues", 
               array3d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleReferenceShapeFunctionGradientValues" );

DECLARE_FIELD( particleCohesiveReferenceSurfaceNormal, 
               "particleCohesiveReferenceSurfaceNormal", 
               array2d< real64 >,
               0.0,
               LEVEL_1,
               WRITE_AND_READ,
               "ParticleCohesiveReferenceSurfaceNormal" );

DECLARE_FIELD( particleCohesiveFieldMapping, 
               "particleCohesiveFieldMapping", 
               array2d< int >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "particleCohesiveFieldMapping" );     

DECLARE_FIELD( particleSubdivideFlag, 
               "particleSubdivideFlag", 
               array1d< int >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "particleSubdivideFlag" );   

DECLARE_FIELD( particleCopyFlag, 
               "particleCopyFlag", 
               array1d< int >,
               -1,
               LEVEL_1,
               WRITE_AND_READ,
               "particleCopyFlag" );

DECLARE_FIELD( particleDomainScaledFlag, 
               "particleDomainScaledFlag", 
               array1d< int >,
               -1,
               LEVEL_1,
               WRITE_AND_READ,
               "particleDomainScaledFlag" );            
}

}

}

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_MPMSOLVERBASEFIELDS_HPP_
