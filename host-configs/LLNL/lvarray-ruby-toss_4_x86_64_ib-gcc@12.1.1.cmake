#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: gcc@=12.1.1
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.3.1/cmake-3.26.3-nz532rvfpaf5lf74zxmplgiobuhol7lu/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-12.1.1-magic/bin/gcc" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-12.1.1-magic/bin/g++" CACHE PATH "")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG" CACHE STRING "")

set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "")

#--------------------------------------------------------------------------------
# Cuda
#--------------------------------------------------------------------------------

set(ENABLE_CUDA OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/chai-git.df7741f1dbbdc5fff5f7d626151fdf1904e62b19_develop-iev43uxnll3jsod7jloak5iwmox6owlv" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/raja-git.4d7fcba55ebc7cb972b7cc9f6778b48e43792ea1_develop-ikig7gipi7bsuxopxfac7uvrup3kd4rx" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/umpire-git.abd729f40064175e999a83d11d6b073dac4c01d2_develop-4t7pxjbfakzclltrwpajtryfkntshmwu" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/camp-git.0f07de4240c42e0b38a8d872a20440cb4b33d9f5_main-fv2sc4dnl7tma5k3wpu5r2tgcinfpz7u" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/caliper-2.11.0-nwqfzu56byuwrrc5bzymfli25niqhdr4" CACHE PATH "")

set(adiak_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/adiak-0.4.0-dimwte7ij4naho6hhcqowidlxnespn7e/lib/cmake/adiak" CACHE PATH "")

#--------------------------------------------------------------------------------
# Documentation
#--------------------------------------------------------------------------------

set(SPHINX_EXECUTABLE "/usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/bin/sphinx-build" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/doxygen-1.8.20-hbxmvlkrwmpt5mvibhths6cdo5rlor3s/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

