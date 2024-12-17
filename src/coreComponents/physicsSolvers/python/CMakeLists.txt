
if( ENABLE_PYGEOSX )
  set( physicsSolvers_headers
       ${physicsSolvers_headers}
       python/PySolver.hpp
       python/PySolverType.hpp
       PARENT_SCOPE )
  set( physicsSolvers_sources
       ${physicsSolvers_sources}
       python/PySolver.cpp
       PARENT_SCOPE )

  set( dependencyList ${dependencyList} Python3::Python pylvarray PARENT_SCOPE)
endif()

