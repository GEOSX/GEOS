# Specify solver headers
set( physicsSolvers_headers
      ${physicsSolvers_headers}
     simplePDE/LaplaceBaseH1.hpp
     simplePDE/LaplaceFEM.hpp
     simplePDE/LaplaceFEMKernels.hpp
     simplePDE/PhaseFieldDamageFEM.hpp
     simplePDE/PhaseFieldDamageFEMKernels.hpp 
     simplePDE/PhaseFieldPressurizedDamageFEMKernels.hpp
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
      ${physicsSolvers_sources}
     simplePDE/LaplaceBaseH1.cpp
     simplePDE/LaplaceFEM.cpp
     simplePDE/PhaseFieldDamageFEM.cpp 
     PARENT_SCOPE )
