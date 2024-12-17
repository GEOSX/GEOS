include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/ruby-clang-14.cmake)

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6-magic CACHE PATH "")
# This is here to note the required flags for using valgrind. These will have to be propagated to the TPL's
#set( CMAKE_CXX_FLAGS "-march=x86-64-v2 -mno-avx512f" CACHE STRING "" FORCE)
# This is what is required for using address sanitizer. This should be put into GeosxOptions.cmake when it is incorporated into the build options.
#set( CMAKE_CXX_FLAGS "-g -O2 -fno-omit-frame-pointer -fsanitize=address" CACHE STRING "" FORCE)

include(${CMAKE_CURRENT_LIST_DIR}/llnl-cpu-base.cmake)
