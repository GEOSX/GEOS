# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     inducedSeismicity/inducedSeismicityFields.hpp
     inducedSeismicity/rateAndStateFields.hpp
     inducedSeismicity/QuasiDynamicEQ.hpp
     inducedSeismicity/SeismicityRate.hpp
     inducedSeismicity/kernels/RateAndStateKernels.hpp
     inducedSeismicity/kernels/SeismicityRateKernels.hpp 
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     inducedSeismicity/QuasiDynamicEQ.cpp 
     inducedSeismicity/SeismicityRate.cpp 
     PARENT_SCOPE )
