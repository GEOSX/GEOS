.. _Dependencies:

Third-party dependencies
========================

GEOS makes use of multiple third-party libraries (TPLs) and tools, some of which are mandatory and some optional.
We only test against specific versions, and sometimes even require development snapshots (specific git commits).
Not all of these guarantee backwards compatibility, so we strongly recommend building with these specific versions.

List of third-party libraries and tools
---------------------------------------

The two tables below lists the dependencies with their specific versions and relevant CMake variables.
Some of these libraries may have their own system prerequisites.

Libraries
~~~~~~~~~

The following libraries are linked to by GEOS:

============= ========== =========================== ============================= =====================================
Name          Version    Enable option               Path variable                 Description
============= ========== =========================== ============================= =====================================
Adiak_        0.2.0      :code:`ENABLE_CALIPER`      :code:`ADIAK_DIR`             Library for collecting metadata from HPC application runs, and distributing that metadata to subscriber tools.
Caliper_      2.4.0      :code:`ENABLE_CALIPER`      :code:`CALIPER_DIR`           Instrumentation and performance profiling library.
conduit_      0.5.0      *mandatory*                 :code:`CONDUIT_DIR`           Simplified Data Exchange for HPC Simulations.
CHAI_         2.2.2      *mandatory*                 :code:`CHAI_DIR`              Copy-hiding array abstraction to automatically migrate data between memory spaces.
RAJA_         0.12.1     *mandatory*                 :code:`RAJA_DIR`              Collection of C++ software abstractions that enable architecture portability for HPC applications.
hdf5_         1.10.5     *mandatory*                 :code:`HDF5_DIR`              High-performance data management and storage suite.
mathpresso_   2015-12-15 :code:`ENABLE_MATHPRESSO`   :code:`MATHPRESSO_DIR`        Mathematical Expression Parser and JIT Compiler.
pugixml_      1.8.0      *mandatory*                 :code:`PUGIXML_DIR`           Light-weight, simple and fast XML parser for C++ with XPath support.
parmetis_     4.0.3      *mandatory* (with MPI)      :code:`PARMETIS_DIR`          Parallel Graph Partitioning library. Should be built with 64-bit :code:`idx_t` type.
suitesparse_  5.8.1      :code:`ENABLE_SUITESPARSE`  :code:`SUITESPARSE_DIR`       A suite of sparse matrix software.
superlu_dist_ 0f6efc3    :code:`ENABLE_SUPERLU_DIST` :code:`SUPERLU_DIST_DIR`      General purpose library for the direct solution of large, sparse, nonsymmetric systems of linear equations.
hypre_        2186a8f    :code:`ENABLE_HYPRE`        :code:`HYPRE_DIR`             Library of high performance preconditioners and solvers for large sparse linear systems on massively parallel computers.
PETSc_        3.13.0     :code:`ENABLE_PETSC`        :code:`PETSC_DIR`             Suite of data structures and routines for the scalable (parallel) solution of scientific applications.
Trilinos_     12.18.1    :code:`ENABLE_TRILINOS`     :code:`TRILINOS_DIR`          Collection of reusable scientific software libraries, known in particular for linear solvers.
silo_         4.10.3     *mandatory*                 :code:`SILO_DIR`              A Mesh and Field I/O Library and Scientific Database.
VTK_          9.0.0-rc3  :code:`ENABLE_VTK`          :code:`VTK_DIR`               Open source software for manipulating and displaying scientific data.
============= ========== =========================== ============================= =====================================

Tools
~~~~~~~~~

The following tools are used as part of the build process to support GEOS development:

============= ========== =========================== ============================= =====================================
Name          Version    Enable option               Path variable                 Description
============= ========== =========================== ============================= =====================================
doxygen_      1.8.20     :code:`ENABLE_DOXYGEN`      :code:`DOXYGEN_EXECUTABLE`    De facto standard tool for generating documentation from annotated C++ sources.
sphinx_       1.8.5      :code:`ENABLE_SPHINX`       :code:`SPHINX_EXECUTABLE`     A tool that makes it easy to create intelligent and beautiful documentation.
uncrustify_   401a409    :code:`ENABLE_UNCRUSTIFY`   :code:`UNCRUSTIFY_EXECUTABLE` A source code beautifier for C, C++, C#, ObjectiveC, D, Java, Pawn and VALA.
============= ========== =========================== ============================= =====================================

.. _Adiak : https://github.com/LLNL/Adiak
.. _Caliper: https://github.com/LLNL/Caliper
.. _conduit: https://github.com/LLNL/conduit
.. _CHAI : https://github.com/LLNL/CHAI
.. _RAJA : https://github.com/LLNL/RAJA
.. _hdf5 : https://portal.hdfgroup.org/display/HDF5/HDF5
.. _mathpresso : https://github.com/kobalicek/mathpresso
.. _pugixml : https://pugixml.org
.. _parmetis : http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview
.. _silo : https://wci.llnl.gov/simulation/computer-codes/silo
.. _VTK : https://vtk.org/
.. _suitesparse : https://people.engr.tamu.edu/davis/suitesparse.html
.. _superlu_dist : https://portal.nersc.gov/project/sparse/superlu
.. _hypre : https://github.com/hypre-space/hypre
.. _PETSc : https://www.mcs.anl.gov/petsc
.. _Trilinos : https://trilinos.github.io
.. _doxygen : https://www.doxygen.nl/index.html
.. _sphinx : https://www.sphinx-doc.org/en/master/
.. _uncrustify : http://uncrustify.sourceforge.net
.. _GoogleTest : https://github.com/google/googletest
.. _GoogleBenchmark : https://github.com/google/benchmark
.. _BLT : https://github.com/LLNL/blt

Some other dependencies (GoogleTest_, GoogleBenchmark_) are provided through BLT_ build system which is embedded in GEOS source.
No actions are needed to build them.

If you would like to create a Docker image with all dependencies, take a look at
`Dockerfiles <https://github.com/GEOS-DEV/thirdPartyLibs/tree/master/docker>`_
that are used in our CI process.

Building bundled dependencies
-----------------------------

To simplify the process of building TPLs, we provide a git repository `thirdPartyLibs <https://github.com/GEOS-DEV/thirdPartyLibs>`_.
It contains source copies of exact TPL versions required and is updated periodically.
It also contains a CMake script for building all TPLs in a single command.

The recommended steps to build TPLs are:

- Create a host-config file that sets all system-specific CMake variables (compiler and library paths, configuration flags, etc.)
  Take a look at `host-config examples <https://github.com/GEOS-DEV/GEOS/blob/develop/host-configs>`_.
- Configure via ``config-build.py`` script:

  .. code-block:: console

     cd thirdPartyLibs
     python scripts/config-build.py --hostconfig=/path/to/host-config.cmake --buildtype=Release --installpath=/path/to/install/dir -DNUM_PROC=8

  where

  * ``--buildpath`` or ``-bp`` is the build directory (by default, created under current).
  * ``--installpath`` or ``-ip`` is the installation directory(wraps ``CMAKE_INSTALL_PREFIX``).
  * ``--buildtype`` or ``-bt`` is a wrapper to the ``CMAKE_BUILD_TYPE`` option.
  * ``--hostconfig`` or ``-hc`` is a path to host-config file.
  * all other command-line options are passed to CMake.

- Run the build:

  .. code-block:: console

     cd <buildpath>
     make

  .. warning::
     Do not provide ``-j`` argument to ``make`` here, since the top-level make only launches sub-project builds.
     Instead use ``-DNUM_PROC`` option above, which is passed to each sub-project's ``make`` command.

You may also run the CMake configure step manually instead of relying on ``config-build.py``.
The full TPL build may take anywhere between 15 minutes and 2 hours, depending on your machine, number of threads and libraries enabled.

.. note::
   An exception from the above pattern, ``sphinx`` is currently not a part of the TPL bundle and must be installed with your Python or package manager.

.. note::
   PETSc build currently downloads `pt-scotch <https://www.labri.fr/perso/pelegrin/scotch/scotch_en.html>`_ from the internet.
   If you do not have access to internet, modify the `./configure` step of petsc in `CMakeLists.txt` and change the ``--download-ptscotch`` option accordingly.
   `pt-scotch` also relies on `bison` and `flex`.

Installing dependencies individually
------------------------------------

You may also install each individual TPL separately, either manually or through a package manager.
This is a more difficult route, since you are responsible for configuring dependencies in a compatible manner.
Again, we strongly recommend using the exact versions listed above, to avoid possible build problems.

You may look at `our TPL CMake script <https://github.com/GEOS-DEV/thirdPartyLibs/blob/master/CMakeLists.txt>`_ to see how we configure TPL builds.
