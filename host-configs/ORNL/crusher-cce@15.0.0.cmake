set(CCE_VERSION 15.0.0)

include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/ORNL/crusher-cce@${CCE_VERSION}.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/crusher-base.cmake)

set( CONDUIT_DIR "${GEOSX_TPL_DIR}/conduit" CACHE PATH "" )
set( HDF5_DIR "${GEOSX_TPL_DIR}/hdf5" CACHE PATH "" )

set(ENABLE_SILO TRUE CACHE BOOL "" )
set(SILO_DIR "${GEOSX_TPL_DIR}/silo" CACHE PATH "" )

# set(ENABLE_VTK TRUE CACHE BOOL "" )
# set(VTK_DIR "${GEOSX_TPL_DIR}/vtk" CACHE PATH "" )
# # these are needed due to how vtk searches for dependencies.. hahahahahaha
# set(GLEW_ROOT "${GEOSX_TPL_DIR}/glew" CACHE PATH "" )
# set(HDF5_ROOT "${GEOSX_TPL_DIR}/hdf5" CACHE PATH "" )
# set(NetCDF_ROOT "${GEOSX_TPL_DIR}/netcdf-c" CACHE PATH "" )
# set(THEORA_ROOT "${GEOSX_TPL_DIR}/libtheora" CACHE PATH "" )
# set(OGG_ROOT "${GEOSX_TPL_DIR}/libogg" CACHE PATH "" )
# set(JsonCpp_ROOT "${GEOSX_TPL_DIR}/jsoncpp" CACHE PATH "" )
# set(GL2PS_ROOT "${GEOSX_TPL_DIR}/gl2ps" CACHE PATH "" )
# set(PNG_ROOT "${GEOSX_TPL_DIR}/libpng" CACHE PATH "" )
# set(pugixml_DIR "${GEOSX_TPL_DIR}/pugixml/lib64/cmake/pugixml" CACHE PATH "" )
# set(SQLite3_ROOT "${GEOSX_TPL_DIR}/sqlite" CACHE PATH "" )
# set(LibPROJ_ROOT "${GEOSX_TPL_DIR}/proj" CACHE PATH "" )
# set(X11_ROOT "${GEOSX_TPL_DIR}/libx11" CACHE PATH "" )
# set(EXPAT_ROOT "${GEOSX_TPL_DIR}/expat" CACHE PATH "" )
# set(double-conversion_ROOT "${GEOSX_TPL_DIR}/double-conversion" CACHE PATH "" )
# set(LZ4_ROOT "${GEOSX_TPL_DIR}/lz4" CACHE PATH "" )
# set(utf8cpp_ROOT "${GEOSX_TPL_DIR}/utf8cpp" CACHE PATH "" )
# set(Eigen3_ROOT "${GEOSX_TPL_DIR}/eigen" CACHE PATH "" )
# set(JPEG_ROOT "${GEOSX_TPL_DIR}/libjpeg-turbo" CACHE PATH "" )
# set(TIFF_ROOT "${GEOSX_TPL_DIR}/libtiff" CACHE PATH "" )
# set(Freetype_ROOT "${GEOSX_TPL_DIR}/freetype" CACHE PATH "" )

set( BLAS_DIR "/opt/rocm-${ROCM_VERSION}/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOSX_TPL_DIR}/pugixml" CACHE PATH "" )
set( FMT_DIR "${GEOSX_TPL_DIR}/fmt" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOSX_TPL_DIR}/suite-sparse" CACHE PATH "" )

# HYPRE options
set( ENABLE_HYPRE_DEVICE "CPU" CACHE STRING "" )
set( ENABLE_HYPRE_MIXINT FALSE CACHE STRING "" )
set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-18-ga592bbd12_cce-15.0.0_rel" CACHE PATH "" ) # CPU SERIAL (WORKS)

# set( ENABLE_HYPRE_MIXINT TRUE CACHE STRING "" )
# set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-200-ge907ce401_cce-15.0.0_mixint_rel/" CACHE PATH "" ) # CPU SERIAL MIXED INT (doesn't work on our side)

# set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-18-ga592bbd12_cce-15.0.0_omp_rel" CACHE PATH "" ) # CPU OPENMP (had openmp runtime issue)

# set( ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "" )
# set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-200-ge907ce401_cce-15.0.0_rocm-5.4.0_umpire-2022.3.0_rel" CACHE PATH "" ) # patched hip version int32
# set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-200-ge907ce401_cce-15.0.0_rocm-5.4.0_mixint_umpire-2022.3.0_rel" CACHE PATH "" ) # patched hip verison int32/64


set( ENABLE_CALIPER ON CACHE BOOL "" FORCE )
set( CALIPER_DIR "${GEOSX_TPL_DIR}/caliper" CACHE PATH "" )