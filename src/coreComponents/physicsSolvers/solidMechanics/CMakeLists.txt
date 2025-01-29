# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     solidMechanics/LogLevelsInfo.hpp
     solidMechanics/SolidMechanicsFields.hpp
     solidMechanics/SolidMechanicsLagrangianFEM.hpp
     solidMechanics/SolidMechanicsLagrangianFEM.hpp
     solidMechanics/SolidMechanicsLagrangianSSLE.hpp
     solidMechanics/kernels/SolidMechanicsLagrangianFEMKernels.hpp
     solidMechanics/SolidMechanicsMPM.hpp
     solidMechanics/MPMSolverFields.hpp
     solidMechanics/kernels/ExplicitFiniteStrain.hpp
     solidMechanics/kernels/ExplicitFiniteStrain_impl.hpp
     solidMechanics/kernels/ExplicitMPM.hpp
     solidMechanics/kernels/ExplicitSmallStrain.hpp
     solidMechanics/kernels/ExplicitSmallStrain_impl.hpp
     solidMechanics/kernels/FixedStressThermoPoromechanics.hpp
     solidMechanics/kernels/FixedStressThermoPoromechanics_impl.hpp
     solidMechanics/kernels/ImplicitSmallStrainNewmark.hpp
     solidMechanics/kernels/ImplicitSmallStrainNewmark_impl.hpp
     solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp
     solidMechanics/kernels/ImplicitSmallStrainQuasiStatic_impl.hpp
     solidMechanics/SolidMechanicsStateReset.hpp
     solidMechanics/SolidMechanicsStatistics.hpp
     solidMechanics/contact/ContactSolverBase.hpp
     solidMechanics/contact/FractureState.hpp
     solidMechanics/contact/ContactFields.hpp
     solidMechanics/contact/SolidMechanicsEmbeddedFractures.hpp
     solidMechanics/contact/SolidMechanicsLagrangeContact.hpp
     solidMechanics/contact/SolidMechanicsLagrangeContactBubbleStab.hpp
     solidMechanics/contact/SolidMechanicsAugmentedLagrangianContact.hpp
     solidMechanics/contact/kernels/SolidMechanicsConformingContactKernelsBase.hpp
     solidMechanics/contact/kernels/SolidMechanicsDisplacementJumpUpdateKernels.hpp
     solidMechanics/contact/kernels/SolidMechanicsEFEMKernelsBase.hpp
     solidMechanics/contact/kernels/SolidMechanicsEFEMKernels.hpp
     solidMechanics/contact/kernels/SolidMechanicsEFEMStaticCondensationKernels.hpp
     solidMechanics/contact/kernels/SolidMechanicsEFEMKernelsHelper.hpp
     solidMechanics/contact/kernels/SolidMechanicsALMKernelsBase.hpp
     solidMechanics/contact/kernels/SolidMechanicsALMKernels.hpp
     solidMechanics/contact/kernels/SolidMechanicsConformingContactKernelsHelper.hpp
     solidMechanics/contact/kernels/SolidMechanicsContactFaceBubbleKernels.hpp
     solidMechanics/contact/kernels/SolidMechanicsLagrangeContactKernels.hpp 
     solidMechanics/contact/LogLevelsInfo.hpp
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     solidMechanics/SolidMechanicsLagrangianFEM.cpp
     solidMechanics/SolidMechanicsLagrangianSSLE.cpp
     solidMechanics/SolidMechanicsMPM.cpp
     solidMechanics/SolidMechanicsStateReset.cpp
     solidMechanics/SolidMechanicsStatistics.cpp
     solidMechanics/contact/ContactSolverBase.cpp
     solidMechanics/contact/SolidMechanicsEmbeddedFractures.cpp
     solidMechanics/contact/SolidMechanicsLagrangeContact.cpp
     solidMechanics/contact/SolidMechanicsLagrangeContactBubbleStab.cpp
     solidMechanics/contact/SolidMechanicsAugmentedLagrangianContact.cpp
     PARENT_SCOPE )

#include( solidMechanics/kernels/SolidMechanicsKernels.cmake)
