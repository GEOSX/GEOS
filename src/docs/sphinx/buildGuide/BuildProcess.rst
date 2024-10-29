.. _BuildProcess:

Building GEOS
==============

Build steps
---------------------

- Create a host-config file that sets all system-specific CMake variables.
  Take a look at `host-config examples <https://github.com/GEOS-DEV/GEOS/blob/develop/host-configs>`_.
  We recommend the same host-config is used for both TPL and GEOS builds.
  In particular, certain options (such as ``ENABLE_MPI`` or ``ENABLE_CUDA``) need to match between the two.

- Provide paths to all enabled TPLs.
  This can be done in one of two ways:

  * Provide each path via a separate CMake variable (see :ref:`Dependencies` for path variable names).
  * If you built TPLs from the ``tplMirror`` repository, you can set ``GEOSX_TPL_DIR`` variable in your host-config to point to the TPL installation path, and

    .. code-block:: cmake

       include("/path/to/GEOS/host-configs/tpls.cmake")

    which will set all the individual TPL paths for you.

- Configure via ``config-build.py`` script:

  .. code-block:: console

     cd GEOS
     python scripts/config-build.py --hostconfig=/path/to/host-config.cmake --buildtype=Release --installpath=/path/to/install/dir

  where

  * ``--buildpath`` or ``-bp`` is the build directory (by default, created under current working dir).
  * ``--installpath`` or ``-ip`` is the installation directory(wraps ``CMAKE_INSTALL_PREFIX``).
  * ``--buildtype`` or ``-bt`` is a wrapper to the ``CMAKE_BUILD_TYPE`` option.
  * ``--hostconfig`` or ``-hc`` is a path to host-config file.
  * all unrecognized options are passed to CMake.

  If ``--buildpath`` is not used, build directory is automatically named ``build-<config-filename-without-extension>-<buildtype>``.
  It is possible to keep automatic naming and change the build root directory with ``--buildrootdir``.
  In that case, build path will be set to ``<buildrootdir>/<config-filename-without-extension>-<buildtype>``.
  Both ``--buildpath`` and ``--buildrootdir`` are incompatible and cannot be used in the same time.
  Same pattern is applicable to install path, with ``--installpath`` and ``--installrootdir`` options.

- Run the build:

  .. code-block:: console

     cd <buildpath>
     make -j $(nproc)

You may also run the CMake configure step manually instead of relying on ``config-build.py``.
A full build typically takes between 10 and 30 minutes, depending on chosen compilers, options and number of cores.

Configuration options
---------------------

Below is a list of CMake configuration options, in addition to TPL options above.
Some options, when enabled, require additional settings (e.g. ``ENABLE_CUDA``).
Please see `host-config examples <https://github.com/GEOS-DEV/GEOS/blob/develop/host-configs>`_.

=============================== ========= ==============================================================================
Option                          Default   Explanation
=============================== ========= ==============================================================================
``ENABLE_MPI``                  ``ON``    Build with MPI (also applies to TPLs)
``ENABLE_OPENMP``               ``OFF``   Build with OpenMP (also applies to TPLs)
``ENABLE_CUDA``                 ``OFF``   Build with CUDA (also applies to TPLs)
``ENABLE_CUDA_NVTOOLSEXT``      ``OFF``   Enable CUDA NVTX user instrumentation (via GEOS_MARK_SCOPE or GEOS_MARK_FUNCTION macros)
``ENABLE_HIP``                  ``OFF``   Build with HIP/ROCM (also applies to TPLs)
``ENABLE_DOCS``                 ``ON``    Build documentation (Sphinx and Doxygen)
``ENABLE_WARNINGS_AS_ERRORS``   ``ON``    Treat all warnings as errors
``ENABLE_PVTPackage``           ``ON``    Enable PVTPackage library (required for compositional flow runs)
``ENABLE_TOTALVIEW_OUTPUT``     ``OFF``   Enables TotalView debugger custom view of GEOS data structures
``ENABLE_COV``                  ``OFF``   Enables code coverage
``GEOS_ENABLE_TESTS``           ``ON``    Enables unit testing targets
``GEOSX_LA_INTERFACE``          ``Hypre`` Choiсe of Linear Algebra backend (Hypre/Petsc/Trilinos)
``GEOSX_BUILD_OBJ_LIBS``        ``ON``    Use CMake Object Libraries build
``GEOSX_BUILD_SHARED_LIBS``     ``OFF``   Build ``geosx_core`` as a shared library instead of static
``GEOSX_PARALLEL_COMPILE_JOBS``           Max. number of compile jobs (when using Ninja), in addition to ``-j`` flag
``GEOSX_PARALLEL_LINK_JOBS``              Max. number of link jobs (when using Ninja), in addition to ``-j`` flag
``GEOSX_INSTALL_SCHEMA``        ``ON``    Enables schema generation and installation
=============================== ========= ==============================================================================
