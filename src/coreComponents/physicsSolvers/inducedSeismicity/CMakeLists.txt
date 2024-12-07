# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     inducedSeismicity/inducedSeismicityFields.hpp
     inducedSeismicity/rateAndStateFields.hpp
     inducedSeismicity/QuasiDynamicEQ.hpp
     inducedSeismicity/ImplicitQDRateAndState.hpp
     inducedSeismicity/SeismicityRate.hpp
     inducedSeismicity/SpringSlider.hpp
     inducedSeismicity/kernels/RateAndStateKernels.hpp
     inducedSeismicity/kernels/SeismicityRateKernels.hpp
     inducedSeismicity/tractionUpdateWrapper/FaultTractionUpdate.hpp
     inducedSeismicity/tractionUpdateWrapper/FaultTractionUpdateBase.hpp
     inducedSeismicity/tractionUpdateWrapper/OneWayCoupledTractionUpdate.hpp
     inducedSeismicity/tractionUpdateWrapper/TractionUpdateFactory.hpp
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     inducedSeismicity/QuasiDynamicEQ.cpp
     inducedSeismicity/ImplicitQDRateAndState.cpp 
     inducedSeismicity/SeismicityRate.cpp
     inducedSeismicity/SpringSlider.cpp
     inducedSeismicity/tractionUpdateWrapper/OneWayCoupledTractionUpdate.cpp
     PARENT_SCOPE )
