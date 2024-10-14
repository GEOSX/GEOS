# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     inducedSeismicity/inducedSeismicityFields.hpp
     inducedSeismicity/SeismicityRate.hpp
     inducedSeismicity/SeismicityRateKernels.hpp 
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     inducedSeismicity/SeismicityRate.cpp 
     PARENT_SCOPE )
