#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: gcc@=12noAVX
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.3.1/cmake-3.26.3-nz532rvfpaf5lf74zxmplgiobuhol7lu/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-12.1.1-magic/bin/gcc" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-12.1.1-magic/bin/g++" CACHE PATH "")

set(CMAKE_CXX_FLAGS "-march=x86-64-v2 -mno-avx512f" CACHE PATH "")

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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/chai-git.df7741f1dbbdc5fff5f7d626151fdf1904e62b19_develop-wtbhnhf2zdlchyvk3xfgwni2wnjxo3js" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/raja-git.4d7fcba55ebc7cb972b7cc9f6778b48e43792ea1_develop-e2mjbgvxgfyunnduvgkuntfkblvmrvjo" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/umpire-git.abd729f40064175e999a83d11d6b073dac4c01d2_develop-ixn6mth7nl44zyt34hztcdnqit6h6aqy" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/camp-git.0f07de4240c42e0b38a8d872a20440cb4b33d9f5_main-mjbymoc6da5vwtgx7kgfcg6x37ilciyg" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/caliper-2.11.0-27xlbuch7tnujplljk76psera2vct5kk" CACHE PATH "")

set(adiak_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/adiak-0.4.0-kgrknuu6drw34t7svfcen5t5u37i4jnf/lib/cmake/adiak" CACHE PATH "")

#--------------------------------------------------------------------------------
# Documentation
#--------------------------------------------------------------------------------

set(SPHINX_EXECUTABLE "/usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/bin/sphinx-build" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/doxygen-1.8.20-c4zarmc366msdoizvmau2bs7n76ob7vo/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

