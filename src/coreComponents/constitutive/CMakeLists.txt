# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 Total, S.A
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
#
#--------------------------------------------------------------------------------------------------

#[[
Package: constitutive

Contains the implementation of constitutive models for fluid and rock properties.
#]]

#
# Specify all headers
#
set( constitutive_headers
     ConstitutiveBase.hpp
     ConstitutiveManager.hpp
     ConstitutivePassThru.hpp
     ConstitutivePassThruHandler.hpp
     ExponentialRelation.hpp
     NullModel.hpp
     capillaryPressure/BrooksCoreyCapillaryPressure.hpp
     capillaryPressure/CapillaryPressureBase.hpp
     capillaryPressure/CapillaryPressureFields.hpp
     capillaryPressure/JFunctionCapillaryPressure.hpp
     capillaryPressure/TableCapillaryPressure.hpp
     capillaryPressure/TableCapillaryPressureHelpers.hpp
     capillaryPressure/VanGenuchtenCapillaryPressure.hpp
     capillaryPressure/capillaryPressureSelector.hpp
     capillaryPressure/layouts.hpp
     contact/BartonBandis.hpp
     contact/CoulombFriction.hpp
     contact/FrictionSelector.hpp
     contact/FrictionBase.hpp
     contact/FrictionlessContact.hpp
     contact/HydraulicApertureBase.hpp
     contact/HydraulicApertureRelationSelector.hpp
     contact/HydraulicApertureTable.hpp
     contact/RateAndStateFriction.hpp
     diffusion/ConstantDiffusion.hpp
     diffusion/DiffusionBase.hpp
     diffusion/DiffusionFields.hpp
     diffusion/DiffusionSelector.hpp
     dispersion/DispersionBase.hpp
     dispersion/DispersionFields.hpp
     dispersion/DispersionSelector.hpp
     dispersion/LinearIsotropicDispersion.hpp
     fluid/multifluid/Layouts.hpp
     fluid/multifluid/MultiFluidSelector.hpp
     fluid/multifluid/MultiFluidBase.hpp
     fluid/multifluid/MultiFluidConstants.hpp
     fluid/multifluid/MultiFluidUtils.hpp
     fluid/multifluid/MultiFluidFields.hpp
     fluid/multifluid/blackOil/BlackOilFluidBase.hpp
     fluid/multifluid/blackOil/BlackOilFluid.hpp
     fluid/multifluid/blackOil/DeadOilFluid.hpp
     fluid/multifluid/blackOil/PVTOData.hpp
     fluid/multifluid/CO2Brine/CO2BrineFluid.hpp
     fluid/multifluid/CO2Brine/PhaseModel.hpp
     fluid/multifluid/CO2Brine/functions/BrineEnthalpy.hpp
     fluid/multifluid/CO2Brine/functions/CO2Enthalpy.hpp
     fluid/multifluid/CO2Brine/functions/CO2EOSSolver.hpp
     fluid/multifluid/CO2Brine/functions/CO2Solubility.hpp
     fluid/multifluid/CO2Brine/functions/CO2SolubilityDuanSun.hpp
     fluid/multifluid/CO2Brine/functions/CO2SolubilitySpycherPruess.hpp
     fluid/multifluid/CO2Brine/functions/EzrokhiBrineDensity.hpp
     fluid/multifluid/CO2Brine/functions/EzrokhiBrineViscosity.hpp
     fluid/multifluid/CO2Brine/functions/FenghourCO2Viscosity.hpp
     fluid/multifluid/CO2Brine/functions/FlashModelBase.hpp
     fluid/multifluid/CO2Brine/functions/NoOpPVTFunction.hpp
     fluid/multifluid/CO2Brine/functions/PhillipsBrineDensity.hpp
     fluid/multifluid/CO2Brine/functions/PhillipsBrineViscosity.hpp
     fluid/multifluid/CO2Brine/functions/PureWaterProperties.hpp
     fluid/multifluid/CO2Brine/functions/PVTFunctionBase.hpp
     fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp
     fluid/multifluid/CO2Brine/functions/SpanWagnerCO2Density.hpp
     fluid/multifluid/CO2Brine/functions/WaterDensity.hpp
     fluid/multifluid/compositional/functions/CompositionalProperties.hpp
     fluid/multifluid/compositional/functions/CompositionalPropertiesImpl.hpp
     fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp
     fluid/multifluid/compositional/functions/FugacityCalculator.hpp
     fluid/multifluid/compositional/functions/KValueInitialization.hpp
     fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp
     fluid/multifluid/compositional/functions/RachfordRice.hpp
     fluid/multifluid/compositional/functions/StabilityTest.hpp
     fluid/multifluid/compositional/models/ComponentProperties.hpp
     fluid/multifluid/compositional/models/CompositionalDensity.hpp
     fluid/multifluid/compositional/models/ConstantViscosity.hpp
     fluid/multifluid/compositional/models/CriticalVolume.hpp
     fluid/multifluid/compositional/models/EquationOfState.hpp
     fluid/multifluid/compositional/models/FunctionBase.hpp
     fluid/multifluid/compositional/models/ImmiscibleWaterFlashModel.hpp
     fluid/multifluid/compositional/models/ImmiscibleWaterParameters.hpp
     fluid/multifluid/compositional/models/LohrenzBrayClarkViscosity.hpp
     fluid/multifluid/compositional/models/LohrenzBrayClarkViscosityImpl.hpp
     fluid/multifluid/compositional/models/NegativeTwoPhaseFlashModel.hpp
     fluid/multifluid/compositional/models/ModelParameters.hpp
     fluid/multifluid/compositional/models/NullModel.hpp
     fluid/multifluid/compositional/models/PhaseModel.hpp
     fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp
     fluid/multifluid/compositional/CompositionalMultiphaseFluidUpdates.hpp
     fluid/multifluid/reactive/ReactiveBrineFluid.hpp
     fluid/multifluid/reactive/ReactiveMultiFluid.hpp
     fluid/multifluid/reactive/ReactiveMultiFluidFields.hpp
     fluid/multifluid/reactive/chemicalReactions/EquilibriumReactions.hpp
     fluid/multifluid/reactive/chemicalReactions/KineticReactions.hpp
     fluid/multifluid/reactive/chemicalReactions/ReactionsBase.hpp
     fluid/singlefluid/CompressibleSinglePhaseFluid.hpp
     fluid/singlefluid/ParticleFluid.hpp
     fluid/singlefluid/ParticleFluidBase.hpp
     fluid/singlefluid/ParticleFluidSelector.hpp
     fluid/singlefluid/ProppantSlurryFluid.hpp
     fluid/singlefluid/SingleFluidBase.hpp
     fluid/singlefluid/SingleFluidFields.hpp
     fluid/singlefluid/SlurryFluidBase.hpp
     fluid/singlefluid/SlurryFluidFields.hpp
     fluid/singlefluid/SingleFluidSelector.hpp
     fluid/singlefluid/SlurryFluidSelector.hpp
     fluid/singlefluid/ThermalCompressibleSinglePhaseFluid.hpp
     permeability/CarmanKozenyPermeability.hpp
     permeability/ConstantPermeability.hpp
     permeability/ExponentialDecayPermeability.hpp
     permeability/ParallelPlatesPermeability.hpp
     permeability/PermeabilityBase.hpp
     permeability/PermeabilityFields.hpp
     permeability/PressurePermeability.hpp
     permeability/ProppantPermeability.hpp
     permeability/SlipDependentPermeability.hpp
     permeability/WillisRichardsPermeability.hpp
     relativePermeability/BrooksCoreyBakerRelativePermeability.hpp
     relativePermeability/BrooksCoreyStone2RelativePermeability.hpp
     relativePermeability/BrooksCoreyRelativePermeability.hpp
     relativePermeability/RelativePermeabilityBase.hpp
     relativePermeability/RelativePermeabilityFields.hpp
     relativePermeability/RelativePermeabilityInterpolators.hpp
     relativePermeability/TableRelativePermeability.hpp
     relativePermeability/TableRelativePermeabilityHelpers.hpp
     relativePermeability/TableRelativePermeabilityHysteresis.hpp
     relativePermeability/VanGenuchtenBakerRelativePermeability.hpp
     relativePermeability/VanGenuchtenStone2RelativePermeability.hpp
     relativePermeability/layouts.hpp
     relativePermeability/RelativePermeabilitySelector.hpp
     solid/CompressibleSolid.hpp
     solid/ProppantSolid.hpp
     solid/CoupledSolidBase.hpp
     solid/CoupledSolid.hpp
     solid/DuvautLionsSolid.hpp
     solid/Damage.hpp
     solid/DamageSpectral.hpp
     solid/DamageSpectralUtilities.hpp
     solid/DamageVolDev.hpp
     solid/DruckerPrager.hpp
     solid/DruckerPragerExtended.hpp
     solid/ModifiedCamClay.hpp
     solid/DelftEgg.hpp
     solid/ElasticIsotropic.hpp
     solid/ElasticIsotropicPressureDependent.hpp
     solid/ElasticTransverseIsotropic.hpp
     solid/ElasticOrthotropic.hpp
     solid/InvariantDecompositions.hpp
     solid/PerfectlyPlastic.hpp
     solid/PorousSolid.hpp
     solid/PropertyConversions.hpp
     solid/SolidBase.hpp
     solid/SolidUtilities.hpp
     solid/SolidInternalEnergy.hpp
     solid/SolidModelDiscretizationOps.hpp
     solid/SolidModelDiscretizationOpsFullyAnisotroipic.hpp
     solid/SolidModelDiscretizationOpsIsotropic.hpp
     solid/SolidModelDiscretizationOpsTransverseIsotropic.hpp
     solid/SolidModelDiscretizationOpsOrthotropic.hpp
     solid/CeramicDamage.hpp
     solid/porosity/PorosityFields.hpp
     solid/porosity/BiotPorosity.hpp
     solid/porosity/PorosityBase.hpp
     solid/porosity/PorosityFields.hpp
     solid/porosity/PressurePorosity.hpp
     solid/porosity/ProppantPorosity.hpp
     thermalConductivity/MultiPhaseConstantThermalConductivity.hpp
     thermalConductivity/MultiPhaseThermalConductivityBase.hpp
     thermalConductivity/MultiPhaseThermalConductivityFields.hpp
     thermalConductivity/MultiPhaseThermalConductivitySelector.hpp
     thermalConductivity/MultiPhaseVolumeWeightedThermalConductivity.hpp
     thermalConductivity/SinglePhaseThermalConductivity.hpp
     thermalConductivity/SinglePhaseThermalConductivityBase.hpp
     thermalConductivity/SinglePhaseThermalConductivityFields.hpp
     thermalConductivity/SinglePhaseThermalConductivitySelector.hpp
     thermalConductivity/ThermalConductivityFields.hpp
   )
#
# Specify all sources
#
set( constitutive_sources
     ConstitutiveBase.cpp
     ConstitutiveManager.cpp
     NullModel.cpp
     capillaryPressure/BrooksCoreyCapillaryPressure.cpp
     capillaryPressure/CapillaryPressureBase.cpp
     capillaryPressure/JFunctionCapillaryPressure.cpp
     capillaryPressure/TableCapillaryPressure.cpp
     capillaryPressure/TableCapillaryPressureHelpers.cpp
     capillaryPressure/VanGenuchtenCapillaryPressure.cpp
     contact/BartonBandis.cpp
     contact/CoulombFriction.cpp
     contact/FrictionBase.cpp
     contact/FrictionlessContact.cpp
     contact/HydraulicApertureBase.cpp
     contact/HydraulicApertureTable.cpp
     contact/RateAndStateFriction.cpp
     diffusion/ConstantDiffusion.cpp
     diffusion/DiffusionBase.cpp
     dispersion/DispersionBase.cpp
     dispersion/LinearIsotropicDispersion.cpp
     fluid/multifluid/MultiFluidBase.cpp
     fluid/multifluid/blackOil/BlackOilFluidBase.cpp
     fluid/multifluid/blackOil/BlackOilFluid.cpp
     fluid/multifluid/blackOil/DeadOilFluid.cpp
     fluid/multifluid/blackOil/PVTOData.cpp
     fluid/multifluid/CO2Brine/CO2BrineFluid.cpp
     fluid/multifluid/CO2Brine/functions/PhillipsBrineDensity.cpp
     fluid/multifluid/CO2Brine/functions/PhillipsBrineViscosity.cpp
     fluid/multifluid/CO2Brine/functions/EzrokhiBrineDensity.cpp
     fluid/multifluid/CO2Brine/functions/EzrokhiBrineViscosity.cpp
     fluid/multifluid/CO2Brine/functions/CO2Solubility.cpp
     fluid/multifluid/CO2Brine/functions/CO2SolubilityDuanSun.cpp
     fluid/multifluid/CO2Brine/functions/CO2SolubilitySpycherPruess.cpp
     fluid/multifluid/CO2Brine/functions/FenghourCO2Viscosity.cpp
     fluid/multifluid/CO2Brine/functions/SpanWagnerCO2Density.cpp
     fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.cpp
     fluid/multifluid/CO2Brine/functions/BrineEnthalpy.cpp
     fluid/multifluid/CO2Brine/functions/CO2Enthalpy.cpp
     fluid/multifluid/CO2Brine/functions/CO2EOSSolver.cpp
     fluid/multifluid/CO2Brine/functions/PureWaterProperties.cpp
     fluid/multifluid/CO2Brine/functions/WaterDensity.cpp
     fluid/multifluid/compositional/models/CompositionalDensity.cpp
     fluid/multifluid/compositional/models/ConstantViscosity.cpp
     fluid/multifluid/compositional/models/CriticalVolume.cpp
     fluid/multifluid/compositional/models/ImmiscibleWaterFlashModel.cpp
     fluid/multifluid/compositional/models/ImmiscibleWaterParameters.cpp
     fluid/multifluid/compositional/models/LohrenzBrayClarkViscosity.cpp
     fluid/multifluid/compositional/models/NegativeTwoPhaseFlashModel.cpp
     fluid/multifluid/compositional/CompositionalMultiphaseFluid.cpp
     fluid/multifluid/compositional/CompositionalMultiphaseFluidUpdates.cpp
     fluid/multifluid/reactive/ReactiveBrineFluid.cpp
     fluid/multifluid/reactive/ReactiveMultiFluid.cpp
     fluid/multifluid/reactive/chemicalReactions/EquilibriumReactions.cpp
     fluid/multifluid/reactive/chemicalReactions/KineticReactions.cpp
     fluid/multifluid/reactive/chemicalReactions/ReactionsBase.cpp
     fluid/singlefluid/CompressibleSinglePhaseFluid.cpp
     fluid/singlefluid/ParticleFluid.cpp
     fluid/singlefluid/ParticleFluidBase.cpp
     fluid/singlefluid/ProppantSlurryFluid.cpp
     fluid/singlefluid/SingleFluidBase.cpp
     fluid/singlefluid/SlurryFluidBase.cpp
     fluid/singlefluid/ThermalCompressibleSinglePhaseFluid.cpp
     permeability/CarmanKozenyPermeability.cpp
     permeability/ConstantPermeability.cpp
     permeability/ExponentialDecayPermeability.cpp
     permeability/ParallelPlatesPermeability.cpp
     permeability/PermeabilityBase.cpp
     permeability/PressurePermeability.cpp
     permeability/ProppantPermeability.cpp
     permeability/SlipDependentPermeability.cpp
     permeability/WillisRichardsPermeability.cpp
     relativePermeability/BrooksCoreyBakerRelativePermeability.cpp
     relativePermeability/BrooksCoreyStone2RelativePermeability.cpp
     relativePermeability/BrooksCoreyRelativePermeability.cpp
     relativePermeability/RelativePermeabilityBase.cpp
     relativePermeability/TableRelativePermeability.cpp
     relativePermeability/TableRelativePermeabilityHelpers.cpp
     relativePermeability/TableRelativePermeabilityHysteresis.cpp
     relativePermeability/VanGenuchtenBakerRelativePermeability.cpp
     relativePermeability/VanGenuchtenStone2RelativePermeability.cpp
     solid/CompressibleSolid.cpp
     solid/CoupledSolidBase.cpp
     solid/ProppantSolid.cpp
     solid/DuvautLionsSolid.cpp
     solid/Damage.cpp
     solid/DamageSpectral.cpp
     solid/DamageVolDev.cpp
     solid/DruckerPrager.cpp
     solid/DruckerPragerExtended.cpp
     solid/ModifiedCamClay.cpp
     solid/DelftEgg.cpp
     solid/ElasticIsotropic.cpp
     solid/ElasticIsotropicPressureDependent.cpp
     solid/ElasticTransverseIsotropic.cpp
     solid/ElasticOrthotropic.cpp
     solid/PerfectlyPlastic.cpp
     solid/PorousSolid.cpp
     solid/SolidBase.cpp
     solid/SolidInternalEnergy.cpp
     solid/CeramicDamage.cpp
     solid/porosity/BiotPorosity.cpp
     solid/porosity/PorosityBase.cpp
     solid/porosity/PressurePorosity.cpp
     solid/porosity/ProppantPorosity.cpp
     thermalConductivity/MultiPhaseConstantThermalConductivity.cpp
     thermalConductivity/MultiPhaseThermalConductivityBase.cpp
     thermalConductivity/MultiPhaseVolumeWeightedThermalConductivity.cpp
     thermalConductivity/SinglePhaseThermalConductivity.cpp
     thermalConductivity/SinglePhaseThermalConductivityBase.cpp
   )

set( dependencyList ${parallelDeps} functions denseLinearAlgebra )

if( ENABLE_PVTPackage )
    set( constitutive_headers
         ${constitutive_headers}
         fluid/multifluid/compositional/CompositionalMultiphaseFluidPVTPackage.hpp )

    set( constitutive_sources
         ${constitutive_sources}
         fluid/multifluid/compositional/CompositionalMultiphaseFluidPVTPackage.cpp
          )

    add_subdirectory( PVTPackage )

    list( APPEND dependencyList PVTPackage )
endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       constitutive
                 SOURCES    ${constitutive_sources}
                 HEADERS    ${constitutive_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( constitutive PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS constitutive LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif( )
