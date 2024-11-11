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

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/hdf5-1.12.1-ijx73yvfrhxnjls4acfy2sdqtiekvyec" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-jwijeuulb6xfflmsmh4snh6qvmmvenve" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/silo-4.11.1-plkkhrhgyyt5d6jhzaeijdeawle3347p" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/pugixml-1.13-7v6zhre6pi7ibnief5tqhmqxowllqo2w" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/vtk-9.3.1-6rs7wpfvmzw7pdsdl6mifdrg7nhuidqd" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/fmt-10.0.0-4wcnmovbi74zd66vdk6xhjo6hefdbxbp" CACHE PATH "")

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

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/metis-5.1.0-jrpyyvbvi5zzxt3mraydgurmfk7pjj77" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/parmetis-4.0.3-i2lkaqn2onxsvx3rcssafh7zkgv4xgcb" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/scotch-7.0.3-fl625kpmhevzvzmuotssykvvhaqovqk3" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-xobiy6gamjkagb3cpnlolghmk7ddded4" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/suite-sparse-5.10.1-fbm7gbwba72t4ggduvynvnbcr43uxxs4" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/trilinos-15.1.1-er7phth5m66nmjkpchfmdqddrbpchqaw" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/hypre-git.06da35b1a4b1066a093bc0c6c48aee12bee74cd4_2.31.0-git.12-gpbyntxpfkcvhoeruewnpjtxjvman2av" CACHE PATH "")

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

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/doxygen-1.8.20-hbxmvlkrwmpt5mvibhths6cdo5rlor3s/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2024-11-11-24_AVX_Testing/ruby-gcc-12_tpls/gcc-12.1.1/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-bsad7cne3ccgu3munuxms52yxxhxeob5/bin/uncrustify" CACHE PATH "")

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

