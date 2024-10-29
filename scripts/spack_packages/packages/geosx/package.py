# Copyright 2013-2023 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *
import warnings

import socket
import os

from os import environ as env
from os.path import join as pjoin

# Tested specs are located at scripts/spack_configs/<$SYS_TYPE>/spack.yaml (e.g. %clang@10.0.1)

# WARNING: +petsc variant is yet to be tested.


def cmake_cache_entry(name, value, comment=""):
    """Generate a string for a cmake cache variable"""

    return 'set(%s "%s" CACHE PATH "%s")\n\n' % (name, value, comment)


def cmake_cache_list(name, value, comment=""):
    """Generate a list for a cmake cache variable"""

    indent = 5 + len(name)
    join_str = '\n' + ' ' * indent
    return 'set(%s %s CACHE STRING "%s")\n\n' % (name, join_str.join(value), comment)


def cmake_cache_string(name, string, comment=""):
    """Generate a string for a cmake cache variable"""

    return 'set(%s "%s" CACHE STRING "%s")\n\n' % (name, string, comment)


def cmake_cache_option(name, boolean_value, comment=""):
    """Generate a string for a cmake configuration option"""

    value = "ON" if boolean_value else "OFF"
    return 'set(%s %s CACHE BOOL "%s")\n\n' % (name, value, comment)


class Geosx(CMakePackage, CudaPackage):
    """GEOSX simulation framework."""

    homepage = "https://github.com/GEOS-DEV/GEOS"
    git = "https://github.com/GEOS-DEV/GEOS.git"

    # GEOSX needs submodules to build, but not necessary to build dependencies
    version('develop', branch='develop')

    # SPHINX_BEGIN_VARIANTS

    variant('shared', default=True, description='Build Shared Libs.')
    variant('caliper', default=True, description='Build Caliper support.')
    variant('vtk', default=True, description='Build VTK support.')
    variant('fesapi', default=False, description='Build fesapi support.')
    variant('trilinos', default=True, description='Build Trilinos support.')
    variant('hypre', default=True, description='Build HYPRE support.')
    variant('petsc', default=False, description='Build PETSc support.')
    variant('scotch', default=True, description='Build Scotch support.')
    variant('uncrustify', default=True, description='Build Uncrustify support.')
    variant('lai',
            default='hypre',
            description='Linear algebra interface.',
            values=('trilinos', 'hypre', 'petsc'),
            multi=False)
    variant('pygeosx', default=True, description='Enable pygeosx.')

    # SPHINX_END_VARIANTS

    # variant('tests', default=True, description='Build tests')
    # variant('benchmarks', default=False, description='Build benchmarks')
    # variant('examples', default=False, description='Build examples')

    variant('docs', default=False, description='Build docs')
    variant('addr2line', default=True,
            description='Add support for addr2line.')

    # SPHINX_BEGIN_DEPENDS

    depends_on('cmake@3.23:', type='build')

    depends_on('blt')

    #
    # Virtual packages
    #
    depends_on('mpi')
    depends_on('blas')
    depends_on('lapack')

    #
    # Performance portability
    #
    depends_on('raja +openmp~examples~exercises~shared')

    depends_on('umpire +c+openmp~examples+fortran~device_alloc~shared')

    depends_on('chai@2023.06.0 +raja+openmp~examples~shared')

    depends_on('camp')

    with when('+cuda'):
        for sm_ in CudaPackage.cuda_arch_values:
            depends_on('raja+cuda cuda_arch={0}'.format(sm_), when='cuda_arch={0}'.format(sm_))
            depends_on('umpire+cuda cuda_arch={0}'.format(sm_), when='cuda_arch={0}'.format(sm_))
            depends_on('chai+cuda cuda_arch={0}'.format(sm_), when='cuda_arch={0}'.format(sm_))
            depends_on('camp+cuda cuda_arch={0}'.format(sm_), when='cuda_arch={0}'.format(sm_))

    #
    # IO
    #
    depends_on('hdf5@1.12.1')
    depends_on('silo@4.11~fortran')

    depends_on('conduit@0.8.2~test~fortran~hdf5_compat')

    depends_on('adiak@0.2.2', when='+caliper')
    depends_on('caliper@2.10.0~gotcha~sampler~libunwind~libdw', when='+caliper')

    depends_on('pugixml@1.13')

    depends_on('fmt@10.0.0 cxxstd=14')
    depends_on('vtk@9.2.6', when='+vtk')

    depends_on('fesapi', when='+fesapi')

    #
    # Math
    #
    depends_on('parmetis@4.0.3+int64')

    depends_on('superlu-dist +int64+openmp')

    depends_on('scotch@7.0.3 +mpi +int64', when='+scotch')

    depends_on('suite-sparse@5.10.1+openmp')

    trilinos_build_options = '+openmp'
    trilinos_packages = '+aztec+stratimikos~amesos2~anasazi~belos~ifpack2~muelu~sacado+thyra'
    depends_on('trilinos@13.4.1 ' + trilinos_build_options + trilinos_packages, when='+trilinos')

    depends_on("hypre +superlu-dist+mixedint+mpi+openmp", when='+hypre~cuda')

    depends_on("hypre +cuda+superlu-dist+mixedint+mpi+openmp+umpire+unified-memory cxxflags='-fPIC'", when='+hypre+cuda')
    with when('+cuda'):
        for sm_ in CudaPackage.cuda_arch_values:
            depends_on('hypre+cuda cuda_arch={0}'.format(sm_), when='cuda_arch={0}'.format(sm_))

    depends_on('petsc@3.13.0~hdf5~hypre+int64', when='+petsc')
    depends_on('petsc+ptscotch', when='+petsc+scotch')

    #
    # Python
    #
    depends_on('python')


    #
    # Dev tools
    #
    depends_on('uncrustify', when='+uncrustify')

    #
    # Documentation
    #
    depends_on('doxygen@1.8.20', when='+docs', type='build')
    depends_on('py-sphinx@1.6.3:', when='+docs', type='build')

    # SPHINX_END_DEPENDS

    #
    # Conflicts
    #
    conflicts('~trilinos lai=trilinos', msg='To use Trilinos as the Linear Algebra Interface you must build it.')
    conflicts('~hypre lai=hypre', msg='To use HYPRE as the Linear Algebra Interface you must build it.')
    conflicts('~petsc lai=petsc', msg='To use PETSc as the Linear Algebra Interface you must build it.')

    # Only phases necessary for building dependencies and generate host configs
    phases = ['geos_hostconfig', 'lvarray_hostconfig']
    #phases = ['hostconfig', 'cmake', 'build', 'install']

    @run_after('build')
    @on_package_attributes(run_tests=True)
    def check(self):
        """
        Searches the CMake-generated Makefile for the target ``test``
        and runs it if found.
        """
        with working_dir(self.build_directory):
            ctest('-V', '--force-new-ctest-process', '-j 1')

    @run_after('build')
    def build_docs(self):
        if '+docs' in self.spec:
            with working_dir(self.build_directory):
                make('docs')

    def _get_sys_type(self, spec):
        sys_type = str(spec.architecture)
        # if on llnl systems, we can use the SYS_TYPE
        if "SYS_TYPE" in env:
            sys_type = env["SYS_TYPE"]
        return sys_type

    def _get_host_config_path(self, spec, lvarray=False):
        var = ''
        if '+cuda' in spec:
            var = '-'.join([var, 'cuda'])

        hostname = socket.gethostname().rstrip('1234567890')

        if lvarray:
            hostname = "lvarray-" + hostname

        host_config_path = "%s-%s-%s%s.cmake" % (hostname, self._get_sys_type(spec), (str(spec.compiler)).replace('=',''), var)

        dest_dir = self.stage.source_path
        host_config_path = os.path.abspath(pjoin(dest_dir, host_config_path))
        return host_config_path

    def geos_hostconfig(self, spec, prefix, py_site_pkgs_dir=None):
        """
        This method creates a 'host-config' file that specifies
        all of the options used to configure and build GEOSX.

        Note:
          The `py_site_pkgs_dir` arg exists to allow a package that
          subclasses this package provide a specific site packages
          dir when calling this function. `py_site_pkgs_dir` should
          be an absolute path or `None`.

          This is necessary because the spack `site_packages_dir`
          var will not exist in the base class. For more details
          on this issue see: https://github.com/spack/spack/issues/6261
        """

        #######################
        # Compiler Info
        #######################
        c_compiler = env["SPACK_CC"]
        cpp_compiler = env["SPACK_CXX"]

        #######################################################################
        # By directly fetching the names of the actual compilers we appear
        # to doing something evil here, but this is necessary to create a
        # 'host config' file that works outside of the spack install env.
        #######################################################################

        sys_type = self._get_sys_type(spec)

        ##############################################
        # Find and record what CMake is used
        ##############################################

        cmake_exe = spec['cmake'].command.path
        cmake_exe = os.path.realpath(cmake_exe)

        host_config_path = self._get_host_config_path(spec)
        with open(host_config_path, "w") as cfg:
            cfg.write("#{0}\n".format("#" * 80))
            cfg.write("# Generated host-config - Edit at own risk!\n")
            cfg.write("#{0}\n".format("#" * 80))

            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# SYS_TYPE: {0}\n".format(sys_type))
            cfg.write("# Compiler Spec: {0}\n".format(spec.compiler))
            cfg.write("# CMake executable path: %s\n" % cmake_exe)
            cfg.write("#{0}\n\n".format("-" * 80))

            #######################
            # Compiler Settings
            #######################

            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# Compilers\n")
            cfg.write("#{0}\n\n".format("-" * 80))
            cfg.write(cmake_cache_entry("CMAKE_C_COMPILER", c_compiler))
            cflags = ' '.join(spec.compiler_flags['cflags'])
            if cflags:
                cfg.write(cmake_cache_entry("CMAKE_C_FLAGS", cflags))

            cfg.write(cmake_cache_entry("CMAKE_CXX_COMPILER", cpp_compiler))
            cxxflags = ' '.join(spec.compiler_flags['cxxflags'])
            if cxxflags:
                cfg.write(cmake_cache_entry("CMAKE_CXX_FLAGS", cxxflags))

            release_flags = "-O3 -DNDEBUG"
            cfg.write(cmake_cache_string("CMAKE_CXX_FLAGS_RELEASE", release_flags))
            reldebinf_flags = "-O3 -g -DNDEBUG"
            cfg.write(cmake_cache_string("CMAKE_CXX_FLAGS_RELWITHDEBINFO", reldebinf_flags))
            debug_flags = "-O0 -g"
            cfg.write(cmake_cache_string("CMAKE_CXX_FLAGS_DEBUG", debug_flags))

            if "%clang arch=linux-rhel7-ppc64le" in spec:
                cfg.write(cmake_cache_entry("CMAKE_EXE_LINKER_FLAGS", "-Wl,--no-toc-optimize"))

            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# CMake Standard\n")
            cfg.write("#{0}\n\n".format("-" * 80))

            cfg.write(cmake_cache_string("BLT_CXX_STD", "c++17"))

            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# MPI\n")
            cfg.write("#{0}\n\n".format("-" * 80))

            cfg.write(cmake_cache_option('ENABLE_MPI', True))
            cfg.write(cmake_cache_entry('MPI_C_COMPILER', spec['mpi'].mpicc))
            cfg.write(cmake_cache_entry('MPI_CXX_COMPILER', spec['mpi'].mpicxx))

            if sys_type in ('linux-rhel7-ppc64le', 'linux-rhel8-ppc64le', 'blueos_3_ppc64le_ib_p9'):
                cfg.write(cmake_cache_option('ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC', True))
                if socket.gethostname().rstrip('1234567890') == "lassen":
                    cfg.write(cmake_cache_entry('MPIEXEC', 'lrun'))
                    cfg.write(cmake_cache_entry('MPIEXEC_NUMPROC_FLAG', '-n'))
                else:
                    cfg.write(cmake_cache_entry('MPIEXEC', 'jsrun'))
                    cfg.write(cmake_cache_list('MPIEXEC_NUMPROC_FLAG', ['-g1', '--bind', 'rs', '-n']))

            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# OpenMP\n")
            cfg.write("#{0}\n\n".format("-" * 80))

            cfg.write(cmake_cache_option('ENABLE_OPENMP', True))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Cuda\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if '+cuda' in spec:
                cfg.write(cmake_cache_option('ENABLE_CUDA', True))
                cfg.write(cmake_cache_entry('CMAKE_CUDA_STANDARD', 17))

                cudatoolkitdir = spec['cuda'].prefix
                cfg.write(cmake_cache_entry('CUDA_TOOLKIT_ROOT_DIR', cudatoolkitdir))
                cudacompiler = '${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc'
                cfg.write(cmake_cache_entry('CMAKE_CUDA_COMPILER', cudacompiler))

                cmake_cuda_flags = ('-restrict --expt-extended-lambda -Werror '
                                    'cross-execution-space-call,reorder,'
                                    'deprecated-declarations')

                archSpecifiers = ('-mtune', '-mcpu', '-march', '-qtune', '-qarch')
                for archSpecifier in archSpecifiers:
                    for compilerArg in spec.compiler_flags['cxxflags']:
                        if compilerArg.startswith(archSpecifier):
                            cmake_cuda_flags += ' -Xcompiler ' + compilerArg

                if not spec.satisfies('cuda_arch=none'):
                    cuda_arch = spec.variants['cuda_arch'].value
                    cmake_cuda_flags += ' -arch sm_{0}'.format(cuda_arch[0])
                    cfg.write(cmake_cache_string('CMAKE_CUDA_ARCHITECTURES', cuda_arch[0]))

                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS', cmake_cuda_flags))

                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS_RELEASE', '-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3 -Xcompiler -mcpu=powerpc64le -Xcompiler -mtune=powerpc64le'))
                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS_RELWITHDEBINFO', '-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}'))
                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS_DEBUG', '-g -G -O0 -Xcompiler -O0'))

            else:
                cfg.write(cmake_cache_option('ENABLE_CUDA', False))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Performance Portability TPLs\n')
            cfg.write('#{0}\n\n'.format('-' * 80))

            cfg.write(cmake_cache_option('ENABLE_CHAI', True))
            cfg.write(cmake_cache_entry('CHAI_DIR', spec['chai'].prefix))

            cfg.write(cmake_cache_entry('RAJA_DIR', spec['raja'].prefix))

            cfg.write(cmake_cache_option('ENABLE_UMPIRE', True))
            cfg.write(cmake_cache_entry('UMPIRE_DIR', spec['umpire'].prefix))

            cfg.write(cmake_cache_entry('CAMP_DIR', spec['camp'].prefix))

            # yapf: disable
            io_tpls = (
                ('hdf5', 'HDF5', True),
                ('conduit', 'CONDUIT', True),
                ('silo', 'SILO', True),
                ('pugixml', 'PUGIXML', True),
                ('vtk', 'VTK', '+vtk' in spec),
                ('fesapi', 'FESAPI', '+fesapi' in spec),
                ('fmt', 'FMT', True)
            )
            # yapf: enable

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# IO TPLs\n')
            cfg.write('#{0}\n\n'.format('-' * 80))

            if '+caliper' in spec:
                cfg.write(cmake_cache_option('ENABLE_CALIPER', True))
                cfg.write(cmake_cache_entry('CALIPER_DIR', spec['caliper'].prefix))
                cfg.write(cmake_cache_entry('adiak_DIR', spec['adiak'].prefix + '/lib/cmake/adiak'))

            for tpl, cmake_name, enable in io_tpls:
                if enable:
                    cfg.write(cmake_cache_entry('{}_DIR'.format(cmake_name), spec[tpl].prefix))
                else:
                    cfg.write(cmake_cache_option('ENABLE_{}'.format(cmake_name), False))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# System Math Libraries\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if '+intel-oneapi-mkl' in spec:
                cfg.write(cmake_cache_option('ENABLE_MKL', True))
                cfg.write(cmake_cache_entry('MKL_INCLUDE_DIRS', spec['intel-oneapi-mkl'].prefix.include))
                cfg.write(cmake_cache_list('MKL_LIBRARIES', spec['intel-oneapi-mkl'].libs))
            elif '+mkl' in spec:
                cfg.write(cmake_cache_option('ENABLE_MKL', True))
                cfg.write(cmake_cache_entry('MKL_INCLUDE_DIRS', spec['intel-mkl'].prefix.include))
                cfg.write(cmake_cache_list('MKL_LIBRARIES', spec['intel-mkl'].libs))
            elif '+essl' in spec:
                cfg.write(cmake_cache_option('ENABLE_ESSL', True))
                cfg.write(cmake_cache_entry('ESSL_INCLUDE_DIRS', spec['essl'].prefix.include))
                cfg.write(cmake_cache_list('ESSL_LIBRARIES', spec['essl'].libs + spec['cuda'].libs))

                cfg.write(cmake_cache_option('FORTRAN_MANGLE_NO_UNDERSCORE', True))
            else:
                cfg.write(cmake_cache_list('BLAS_LIBRARIES', spec['blas'].libs))
                cfg.write(cmake_cache_list('LAPACK_LIBRARIES', spec['lapack'].libs))

            # yapf: disable
            math_tpls = (
                ('metis', 'METIS', True),
                ('parmetis', 'PARMETIS', True),
                ('scotch', 'SCOTCH', '+scotch' in spec),
                ('superlu-dist', 'SUPERLU_DIST', True),
                ('suite-sparse', 'SUITESPARSE', True),
                ('trilinos', 'TRILINOS', '+trilinos' in spec),
                ('hypre', 'HYPRE', '+hypre' in spec),
                ('petsc', 'PETSC', '+petsc' in spec)
            )
            # yapf: enable

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Math TPLs\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            for tpl, cmake_name, enable in math_tpls:
                if enable:
                    cfg.write(cmake_cache_entry('{}_DIR'.format(cmake_name), spec[tpl].prefix))

                    if tpl == 'hypre' and '+cuda' in spec:
                        cfg.write(cmake_cache_string('ENABLE_HYPRE_DEVICE', "CUDA"))
                    elif tpl == 'hypre' and '+hip' in spec:
                        cfg.write(cmake_cache_string('ENABLE_HYPRE_DEVICE', "HIP"))
                else:
                    cfg.write(cmake_cache_option('ENABLE_{}'.format(cmake_name), False))

            if '+caliper' in spec and '+hypre' in spec and '+cuda' not in spec:
                cfg.write(cmake_cache_option('ENABLE_CALIPER_HYPRE', True))

            if 'lai=trilinos' in spec:
                cfg.write(cmake_cache_entry('GEOS_LA_INTERFACE', 'Trilinos'))
            if 'lai=hypre' in spec:
                cfg.write(cmake_cache_entry('GEOS_LA_INTERFACE', 'Hypre'))
            if 'lai=petsc' in spec:
                cfg.write(cmake_cache_entry('GEOS_LA_INTERFACE', 'Petsc'))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Python\n')
            cfg.write('#{0}\n\n'.format('-' * 80))

            cfg.write(cmake_cache_entry('Python3_ROOT_DIR', os.path.join(spec['python'].prefix)))
            cfg.write(cmake_cache_entry('Python3_EXECUTABLE', os.path.join(spec['python'].prefix.bin, 'python3')))

            if '+pygeosx' in spec:
                cfg.write(cmake_cache_option('ENABLE_PYGEOSX', True))
            else:
                cfg.write(cmake_cache_option('ENABLE_PYGEOSX', False))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Documentation\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if '+docs' in spec:
                sphinx_bin_dir = spec['py-sphinx'].prefix.bin
                cfg.write(cmake_cache_entry('SPHINX_EXECUTABLE', os.path.join(sphinx_bin_dir, 'sphinx-build')))

                doxygen_bin_dir = spec['doxygen'].prefix.bin
                cfg.write(cmake_cache_entry('DOXYGEN_EXECUTABLE', os.path.join(doxygen_bin_dir, 'doxygen')))
            else:
                cfg.write(cmake_cache_option('ENABLE_DOCS', False))
                cfg.write(cmake_cache_option('ENABLE_DOXYGEN', False))
                cfg.write(cmake_cache_option('ENABLE_SPHINX', False))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Development tools\n')
            cfg.write('#{0}\n\n'.format('-' * 80))

            cfg.write(cmake_cache_option('ENABLE_UNCRUSTIFY', '+uncrustify' in spec))
            if '+uncrustify' in spec:
                cfg.write(
                    cmake_cache_entry('UNCRUSTIFY_EXECUTABLE', os.path.join(spec['uncrustify'].prefix.bin, 'uncrustify')))

            if '+addr2line' in spec:
                cfg.write('#{0}\n'.format('-' * 80))
                cfg.write('# addr2line\n')
                cfg.write('#{0}\n\n'.format('-' * 80))
                cfg.write(cmake_cache_option('ENABLE_ADDR2LINE', True))
                cfg.write(cmake_cache_entry('ADDR2LINE_EXEC ', '/usr/bin/addr2line'))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Other\n')
            cfg.write('#{0}\n\n'.format('-' * 80))

            cfg.write(cmake_cache_option('ENABLE_MATHPRESSO', False))
            cfg.write(cmake_cache_option('ENABLE_XML_UPDATES', False))

            # ATS
            # Lassen
            if sys_type in ('blueos_3_ppc64le_ib_p9'):
                cfg.write(cmake_cache_string('ATS_ARGUMENTS', '--ats jsrun_omp --ats jsrun_bind=packed'))
            # Quartz
            if sys_type in ('toss_4_x86_64_ib'):
                cfg.write(cmake_cache_string('ATS_ARGUMENTS', '--machine slurm36'))

    def lvarray_hostconfig(self, spec, prefix, py_site_pkgs_dir=None):
        """
        This method creates a 'host-config' file that specifies
        all of the options used to configure and build LvArray.

        Note:
          The `py_site_pkgs_dir` arg exists to allow a package that
          subclasses this package provide a specific site packages
          dir when calling this function. `py_site_pkgs_dir` should
          be an absolute path or `None`.

          This is necessary because the spack `site_packages_dir`
          var will not exist in the base class. For more details
          on this issue see: https://github.com/spack/spack/issues/6261
        """

        #######################
        # Compiler Info
        #######################
        c_compiler = env["SPACK_CC"]
        cpp_compiler = env["SPACK_CXX"]

        #######################################################################
        # By directly fetching the names of the actual compilers we appear
        # to doing something evil here, but this is necessary to create a
        # 'host config' file that works outside of the spack install env.
        #######################################################################

        sys_type = self._get_sys_type(spec)

        ##############################################
        # Find and record what CMake is used
        ##############################################

        cmake_exe = spec['cmake'].command.path
        cmake_exe = os.path.realpath(cmake_exe)

        host_config_path = self._get_host_config_path(spec, lvarray=True)
        with open(host_config_path, "w") as cfg:
            cfg.write("#{0}\n".format("#" * 80))
            cfg.write("# Generated host-config - Edit at own risk!\n")
            cfg.write("#{0}\n".format("#" * 80))

            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# SYS_TYPE: {0}\n".format(sys_type))
            cfg.write("# Compiler Spec: {0}\n".format(spec.compiler))
            cfg.write("# CMake executable path: %s\n" % cmake_exe)
            cfg.write("#{0}\n\n".format("-" * 80))

            #######################
            # Compiler Settings
            #######################

            cfg.write("#{0}\n".format("-" * 80))
            cfg.write("# Compilers\n")
            cfg.write("#{0}\n\n".format("-" * 80))
            cfg.write(cmake_cache_entry("CMAKE_C_COMPILER", c_compiler))
            cflags = ' '.join(spec.compiler_flags['cflags'])
            if cflags:
                cfg.write(cmake_cache_entry("CMAKE_C_FLAGS", cflags))

            cfg.write(cmake_cache_entry("CMAKE_CXX_COMPILER", cpp_compiler))
            cxxflags = ' '.join(spec.compiler_flags['cxxflags'])
            if cxxflags:
                cfg.write(cmake_cache_entry("CMAKE_CXX_FLAGS", cxxflags))

            release_flags = "-O3 -DNDEBUG"
            cfg.write(cmake_cache_string("CMAKE_CXX_FLAGS_RELEASE", release_flags))
            reldebinf_flags = "-O3 -g -DNDEBUG"
            cfg.write(cmake_cache_string("CMAKE_CXX_FLAGS_RELWITHDEBINFO", reldebinf_flags))
            debug_flags = "-O0 -g"
            cfg.write(cmake_cache_string("CMAKE_CXX_FLAGS_DEBUG", debug_flags))

            if "%clang arch=linux-rhel7-ppc64le" in spec:
                cfg.write(cmake_cache_entry("CMAKE_EXE_LINKER_FLAGS", "-Wl,--no-toc-optimize"))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Cuda\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if '+cuda' in spec:
                cfg.write(cmake_cache_option('ENABLE_CUDA', True))
                cfg.write(cmake_cache_entry('CMAKE_CUDA_STANDARD', 17))

                cudatoolkitdir = spec['cuda'].prefix
                cfg.write(cmake_cache_entry('CUDA_TOOLKIT_ROOT_DIR', cudatoolkitdir))
                cudacompiler = '${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc'
                cfg.write(cmake_cache_entry('CMAKE_CUDA_COMPILER', cudacompiler))

                cmake_cuda_flags = ('-restrict --expt-extended-lambda -Werror '
                                    'cross-execution-space-call,reorder,'
                                    'deprecated-declarations')

                archSpecifiers = ('-mtune', '-mcpu', '-march', '-qtune', '-qarch')
                for archSpecifier in archSpecifiers:
                    for compilerArg in spec.compiler_flags['cxxflags']:
                        if compilerArg.startswith(archSpecifier):
                            cmake_cuda_flags += ' -Xcompiler ' + compilerArg

                if not spec.satisfies('cuda_arch=none'):
                    cuda_arch = spec.variants['cuda_arch'].value
                    cmake_cuda_flags += ' -arch sm_{0}'.format(cuda_arch[0])
                    cfg.write(cmake_cache_string('CMAKE_CUDA_ARCHITECTURES', cuda_arch[0]))

                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS', cmake_cuda_flags))

                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS_RELEASE', '-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3 -Xcompiler -mcpu=powerpc64le -Xcompiler -mtune=powerpc64le'))
                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS_RELWITHDEBINFO', '-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}'))
                cfg.write(cmake_cache_string('CMAKE_CUDA_FLAGS_DEBUG', '-g -G -O0 -Xcompiler -O0'))

            else:
                cfg.write(cmake_cache_option('ENABLE_CUDA', False))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Performance Portability TPLs\n')
            cfg.write('#{0}\n\n'.format('-' * 80))

            cfg.write(cmake_cache_option('ENABLE_CHAI', True))
            cfg.write(cmake_cache_entry('CHAI_DIR', spec['chai'].prefix))

            cfg.write(cmake_cache_entry('RAJA_DIR', spec['raja'].prefix))

            cfg.write(cmake_cache_option('ENABLE_UMPIRE', True))
            cfg.write(cmake_cache_entry('UMPIRE_DIR', spec['umpire'].prefix))

            cfg.write(cmake_cache_entry('CAMP_DIR', spec['camp'].prefix))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# IO TPLs\n')
            cfg.write('#{0}\n\n'.format('-' * 80))

            if '+caliper' in spec:
                cfg.write(cmake_cache_option('ENABLE_CALIPER', True))
                cfg.write(cmake_cache_entry('CALIPER_DIR', spec['caliper'].prefix))
                cfg.write(cmake_cache_entry('adiak_DIR', spec['adiak'].prefix + '/lib/cmake/adiak'))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Documentation\n')
            cfg.write('#{0}\n\n'.format('-' * 80))
            if '+docs' in spec:
                sphinx_bin_dir = spec['py-sphinx'].prefix.bin
                cfg.write(cmake_cache_entry('SPHINX_EXECUTABLE', os.path.join(sphinx_bin_dir, 'sphinx-build')))

                doxygen_bin_dir = spec['doxygen'].prefix.bin
                cfg.write(cmake_cache_entry('DOXYGEN_EXECUTABLE', os.path.join(doxygen_bin_dir, 'doxygen')))
            else:
                cfg.write(cmake_cache_option('ENABLE_DOXYGEN', False))
                cfg.write(cmake_cache_option('ENABLE_SPHINX', False))

            cfg.write('#{0}\n'.format('-' * 80))
            cfg.write('# Development tools\n')
            cfg.write('#{0}\n\n'.format('-' * 80))

            if '+addr2line' in spec:
                cfg.write('#{0}\n'.format('-' * 80))
                cfg.write('# addr2line\n')
                cfg.write('#{0}\n\n'.format('-' * 80))
                cfg.write(cmake_cache_option('ENABLE_ADDR2LINE', True))

    def cmake_args(self):
        pass
        # spec = self.spec
        # host_config_path = self._get_host_config_path(spec)

        # options = []
        # options.extend(['-C', host_config_path])

        # # Shared libs
        # options.append(self.define_from_variant('BUILD_SHARED_LIBS', 'shared'))

        # if '~tests~examples~benchmarks' in spec:
        #     options.append('-DGEOS_ENABLE_TESTS=OFF')
        # else:
        #     options.append('-DGEOS_ENABLE_TESTS=ON')

        # if '~test' in spec:
        #     options.append('-DDISABLE_UNIT_TESTS=ON')
        # elif "+tests" in spec and ('%intel' in spec or '%xl' in spec):
        #     warnings.warn('The LvArray unit tests take an excessive amount of'
        #                   ' time to build with the Intel or IBM compilers.')

        # options.append(self.define_from_variant('ENABLE_EXAMPLES', 'examples'))
        # options.append(self.define_from_variant('ENABLE_BENCHMARKS',
        #                                         'benchmarks'))
        # options.append(self.define_from_variant('ENABLE_DOCS', 'docs'))

        # return options
