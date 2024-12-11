# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     inducedSeismicity/ExplicitQDRateAndState.hpp
     inducedSeismicity/ImplicitQDRateAndState.hpp
     inducedSeismicity/inducedSeismicityFields.hpp
     inducedSeismicity/QuasiDynamicEarthQuake.hpp
     inducedSeismicity/rateAndStateFields.hpp
     inducedSeismicity/SeismicityRate.hpp
     inducedSeismicity/SpringSlider.hpp
     inducedSeismicity/kernels/ExplicitRateAndStateKernels.hpp
     inducedSeismicity/kernels/ImplicitRateAndStateKernels.hpp
     inducedSeismicity/kernels/RateAndStateKernelsBase.hpp
     inducedSeismicity/kernels/SeismicityRateKernels.hpp
     inducedSeismicity/tractionUpdateWrapper/FaultTractionUpdate.hpp
     inducedSeismicity/tractionUpdateWrapper/FaultTractionUpdateBase.hpp
     inducedSeismicity/tractionUpdateWrapper/OneWayCoupledTractionUpdate.hpp
     inducedSeismicity/tractionUpdateWrapper/TractionUpdateFactory.hpp
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     inducedSeismicity/ExplicitQDRateAndState.cpp 
     inducedSeismicity/ImplicitQDRateAndState.cpp 
     inducedSeismicity/QuasiDynamicEarthQuake.cpp
     inducedSeismicity/SeismicityRate.cpp
     inducedSeismicity/SpringSlider.cpp
     inducedSeismicity/tractionUpdateWrapper/OneWayCoupledTractionUpdate.cpp
     PARENT_SCOPE )
