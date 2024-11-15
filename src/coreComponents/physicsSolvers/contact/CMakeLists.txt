# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     contact/ContactSolverBase.hpp
     contact/ContactFields.hpp
     contact/SolidMechanicsEFEMKernelsBase.hpp
     contact/SolidMechanicsEFEMKernels.hpp
     contact/SolidMechanicsEFEMStaticCondensationKernels.hpp
     contact/SolidMechanicsEFEMKernelsHelper.hpp
     contact/SolidMechanicsEmbeddedFractures.hpp
     contact/SolidMechanicsLagrangeContact.hpp
     contact/SolidMechanicsAugmentedLagrangianContact.hpp
     contact/SolidMechanicsALMKernelsBase.hpp
     contact/SolidMechanicsALMKernels.hpp
     contact/SolidMechanicsALMKernelsHelper.hpp
     contact/SolidMechanicsALMJumpUpdateKernels.hpp
     contact/SolidMechanicsALMBubbleKernels.hpp 
     contact/LogLevelsInfo.hpp
     PARENT_SCOPE )


# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     contact/ContactSolverBase.cpp
     contact/SolidMechanicsEmbeddedFractures.cpp
     contact/SolidMechanicsLagrangeContact.cpp
     contact/SolidMechanicsAugmentedLagrangianContact.cpp
     PARENT_SCOPE )