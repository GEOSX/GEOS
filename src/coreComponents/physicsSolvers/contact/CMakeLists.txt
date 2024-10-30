# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     contact/ContactSolverBase.hpp
     contact/ContactFields.hpp
     contact/SolidMechanicsEmbeddedFractures.hpp
     contact/SolidMechanicsLagrangeContact.hpp
     contact/SolidMechanicsLagrangeContactBubbleStab.hpp
     contact/SolidMechanicsAugmentedLagrangianContact.hpp
     contact/kernels/SolidMechanicsConformingContactKernelsBase.hpp
     contact/kernels/SolidMechanicsDisplacementJumpUpdateKernels.hpp
     contact/kernels/SolidMechanicsEFEMKernelsBase.hpp
     contact/kernels/SolidMechanicsEFEMKernels.hpp
     contact/kernels/SolidMechanicsEFEMStaticCondensationKernels.hpp
     contact/kernels/SolidMechanicsEFEMKernelsHelper.hpp
     contact/kernels/SolidMechanicsALMKernelsBase.hpp
     contact/kernels/SolidMechanicsALMKernels.hpp
     contact/kernels/SolidMechanicsALMKernelsHelper.hpp
     contact/kernels/SolidMechanicsContactFaceBubbleKernels.hpp
     contact/kernels/SolidMechanicsLagrangeContactKernels.hpp 
     contact/LogLevelsInfo.hpp
     PARENT_SCOPE )


# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     contact/ContactSolverBase.cpp
     contact/SolidMechanicsEmbeddedFractures.cpp
     contact/SolidMechanicsLagrangeContact.cpp
     contact/SolidMechanicsLagrangeContactBubbleStab.cpp
     contact/SolidMechanicsAugmentedLagrangianContact.cpp
     PARENT_SCOPE )