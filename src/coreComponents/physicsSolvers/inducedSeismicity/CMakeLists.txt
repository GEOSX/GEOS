# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     inducedSeismicity/ExplicitQDRateAndState.hpp
     inducedSeismicity/ImplicitQDRateAndState.hpp
     inducedSeismicity/inducedSeismicityFields.hpp
     inducedSeismicity/QDRateAndStateBase.hpp
     inducedSeismicity/QuasiDynamicEarthQuake.hpp
     inducedSeismicity/rateAndStateFields.hpp
     inducedSeismicity/SeismicityRate.hpp
     inducedSeismicity/SpringSlider.hpp
     inducedSeismicity/kernels/ExplicitRateAndStateKernels.hpp
     inducedSeismicity/kernels/ImplicitRateAndStateKernels.hpp
     inducedSeismicity/kernels/RateAndStateKernelsBase.hpp
     inducedSeismicity/kernels/SeismicityRateKernels.hpp
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     inducedSeismicity/ExplicitQDRateAndState.cpp 
     inducedSeismicity/ImplicitQDRateAndState.cpp 
     inducedSeismicity/QDRateAndStateBase.cpp
     inducedSeismicity/QuasiDynamicEarthQuake.cpp
     inducedSeismicity/SeismicityRate.cpp
     inducedSeismicity/SpringSlider.cpp
     PARENT_SCOPE )
