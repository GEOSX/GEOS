# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     wavePropagation/shared/WaveSolverBase.hpp
     wavePropagation/shared/WaveSolverUtils.hpp
     wavePropagation/shared/PrecomputeSourcesAndReceiversKernel.hpp
     wavePropagation/sem/acoustic/shared/AcousticFields.hpp
     wavePropagation/sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEM.hpp
     wavePropagation/sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEMKernel.hpp
     wavePropagation/sem/acoustic/shared/AcousticTimeSchemeSEMKernel.hpp
     wavePropagation/sem/acoustic/shared/AcousticMatricesSEMKernel.hpp
     wavePropagation/sem/acoustic/secondOrderEqn/isotropic/AcousticPMLSEMKernel.hpp
     wavePropagation/sem/elastic/shared/ElasticFields.hpp
     wavePropagation/sem/elastic/secondOrderEqn/isotropic/ElasticWaveEquationSEM.hpp
     wavePropagation/sem/elastic/secondOrderEqn/isotropic/ElasticWaveEquationSEMKernel.hpp
     wavePropagation/sem/elastic/shared/ElasticTimeSchemeSEMKernel.hpp
     wavePropagation/sem/elastic/firstOrderEqn/isotropic/ElasticFirstOrderWaveEquationSEM.hpp
     wavePropagation/sem/elastic/firstOrderEqn/isotropic/ElasticFirstOrderWaveEquationSEMKernel.hpp
     wavePropagation/sem/acoustic/firstOrderEqn/isotropic/AcousticFirstOrderWaveEquationSEM.hpp
     wavePropagation/sem/acoustic/firstOrderEqn/isotropic/AcousticFirstOrderWaveEquationSEMKernel.hpp
     wavePropagation/sem/elastic/shared/ElasticMatricesSEMKernel.hpp
     wavePropagation/sem/elastic/secondOrderEqn/anisotropic/ElasticVTIFields.hpp
     wavePropagation/sem/elastic/secondOrderEqn/anisotropic/ElasticVTIWaveEquationSEMKernel.hpp
     wavePropagation/sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIFields.hpp
     wavePropagation/sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIWaveEquationSEM.hpp
     wavePropagation/sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIWaveEquationSEMKernel.hpp
     wavePropagation/sem/acoustoelastic/secondOrderEqn/isotropic/AcoustoElasticFields.hpp
     wavePropagation/sem/acoustoelastic/secondOrderEqn/isotropic/AcoustoElasticTimeSchemeSEMKernel.hpp
     wavePropagation/sem/acoustoelastic/secondOrderEqn/isotropic/AcousticElasticWaveEquationSEM.hpp
     wavePropagation/sem/acoustoelastic/secondOrderEqn/isotropic/AcousticElasticWaveEquationSEMKernel.hpp 
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     wavePropagation/shared/WaveSolverBase.cpp
     wavePropagation/sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEM.cpp
     wavePropagation/sem/elastic/secondOrderEqn/isotropic/ElasticWaveEquationSEM.cpp
     wavePropagation/sem/elastic/firstOrderEqn/isotropic/ElasticFirstOrderWaveEquationSEM.cpp
     wavePropagation/sem/acoustic/firstOrderEqn/isotropic/AcousticFirstOrderWaveEquationSEM.cpp
     wavePropagation/sem/acoustic/secondOrderEqn/anisotropic/AcousticVTIWaveEquationSEM.cpp
     wavePropagation/sem/acoustoelastic/secondOrderEqn/isotropic/AcousticElasticWaveEquationSEM.cpp 
     PARENT_SCOPE)
