# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
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
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     solidMechanics/SolidMechanicsLagrangianFEM.cpp
     solidMechanics/SolidMechanicsLagrangianSSLE.cpp
     solidMechanics/SolidMechanicsMPM.cpp
     solidMechanics/SolidMechanicsStateReset.cpp
     solidMechanics/SolidMechanicsStatistics.cpp
     PARENT_SCOPE )

#include( solidMechanics/kernels/SolidMechanicsKernels.cmake)