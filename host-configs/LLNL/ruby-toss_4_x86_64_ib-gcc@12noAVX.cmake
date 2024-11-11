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
# CMake Standard
#--------------------------------------------------------------------------------

set(BLT_CXX_STD "c++17" CACHE STRING "")

#--------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-12.1.1-magic/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-12.1.1-magic/bin/mpicxx" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

#--------------------------------------------------------------------------------
# OpenMP
#--------------------------------------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

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

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/hdf5-1.12.1-chdwkirwhv5wph3ui6yuenspycyu4o5g" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-dc5zun6rptqzosk3ntitczne2umvfvoe" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/silo-4.11.1-fd4s2c75ynzoy36sqtvogi73fjqloi6m" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/pugixml-1.13-lzydzfiimopvd46asgg6wsswe4tkzbvm" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/vtk-9.3.1-rtmsf25pdwhvyygjqmzsezgoomzqf6fi" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/fmt-10.0.0-na2mo7xln32wmq2hjrmbdyzacr7yqxev" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(ENABLE_MKL ON CACHE BOOL "")

set(MKL_INCLUDE_DIRS "/usr/tce/packages/mkl/mkl-2022.1.0/include" CACHE PATH "")

set(MKL_LIBRARIES /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_gf_lp64.so
                  /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_gnu_thread.so
                  /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_core.so
                  /lib/../lib64/libomp.so
                  /lib64/libpthread.so
                  /lib64/libm.so
                  /lib64/libdl.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/metis-5.1.0-klusyacs6uicwrr6t2uh4inudrz25oog" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/parmetis-4.0.3-v5legva4mfqqce3sdgesl4th77brabgp" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/scotch-7.0.3-o4ndyzzae5dofc4ntxfnw3ujd2qluj2d" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-rxhoos5q7wabnuictdonpxdyu3y263du" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/suite-sparse-5.10.1-webltrnxup4nbupu53no5257splerdr6" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/trilinos-15.1.1-khqsr54db6lh4ol6ynizxuvbum5in3ny" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/hypre-git.06da35b1a4b1066a093bc0c6c48aee12bee74cd4_2.31.0-git.12-hm4ko3o34l57h47gia7zzzccn7qmsytf" CACHE PATH "")

set(ENABLE_PETSC OFF CACHE BOOL "")

set(ENABLE_CALIPER_HYPRE ON CACHE BOOL "")

set(GEOS_LA_INTERFACE "Hypre" CACHE STRING "")

#--------------------------------------------------------------------------------
# Python
#--------------------------------------------------------------------------------

set(Python3_ROOT_DIR "/usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/" CACHE PATH "")

set(Python3_EXECUTABLE "/usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/bin/python3" CACHE PATH "")

set(ENABLE_PYGEOSX ON CACHE BOOL "")

#--------------------------------------------------------------------------------
# Documentation
#--------------------------------------------------------------------------------

set(SPHINX_EXECUTABLE "/usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/bin/sphinx-build" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/doxygen-1.8.20-c4zarmc366msdoizvmau2bs7n76ob7vo/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12noAVX_tpls/gcc-12noAVX/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-2mqopy2akzlsxtxeg27byz3abovta3nh/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

set(ADDR2LINE_EXEC  "/usr/bin/addr2line" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO OFF CACHE BOOL "")

set(ENABLE_XML_UPDATES OFF CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm56" CACHE STRING "")

