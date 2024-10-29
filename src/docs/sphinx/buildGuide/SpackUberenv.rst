.. _SpackUberenv:

Spack and Uberenv
=================

GEOS is transitioning to a new `Spack <https://github.com/spack/spack/>`_ and `Uberenv <https://github.com/LLNL/uberenv>`_ system for building our dependencies. We refer the reader to the `Spack documentation <https://spack.readthedocs.io/en/latest/index.html>`_ and `Uberenv documentation <https://uberenv.readthedocs.io/en/latest/>`_, in particular the Spack documentation for `specs and dependencies <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_, `manual compiler configuration <https://spack.readthedocs.io/en/latest/getting_started.html?highlight=compilers.yaml#manual-compiler-configuration>`_ and `external packages <https://spack.readthedocs.io/en/latest/build_settings.html#external-packages>`_ are worth reading.

Building the dependencies can be as simple as running

.. code-block:: console

    ./scripts/uberenv/uberenv.py

This will create a directory ``uberenv_libs`` (or a directory name you specify by adding ``--prefix directory-name``) in the current working directory, clone Spack into ``uberenv_libs/spack`` and install the dependencies into ``uberenv_libs/system_dependent_path``. It will then spit out a host-config file in the current directory which you can use to build GEOS. While the above command **should** work on every system, it **should never be used**. Invoked as such, Spack will ignore any system libraries you have installed and will go down a rabbit hole building dependencies. Furthermore this does not allow you to choose the compiler to build. Both of these are easily solved by creating a directory with a ``spack.yaml``.

To prevent this from happening you'll need to create a directory with a ``spack.yaml`` file. You can find working examples for commonly used systems in `scripts/spack_configs <https://github.com/GEOS-DEV/GEOS/tree/develop/scripts/spack_configs>`_.

Once you have these files setup you can run Uberenv again and instruct it to use them with. If for instance you added Clang 10.0.1 to the ``spack.yaml`` file the your command would look something like this:

.. code-block:: console

    ./scripts/uberenv/uberenv.py --spack-config-dir=/path/to/your/config/directory/ --spec="%clang@14.0.6"

It is worth noting that GEOS has `two project json files <https://uberenv.readthedocs.io/en/latest/#project-configuration>`_ (``.uberenv_config.json`` and ``scripts/pygeosx_configs/pygeosx.json``) and two configuration directories for LC systems (``scripts/spack_configs`` and ``scripts/pygeosx_configs``). The ``.uberenv_config.json`` project json file and ``scripts/spack_configs`` directory is for building GEOS. The ``scripts/pygeosx_configs/pygeosx.json`` project json file and ``scripts/pygeosx_configs`` directory is for building ``pygeosx``. This is because ``pygeosx`` has a separate list of required compilers and packages to build from GEOS (e.g. ``pygeosx``'s numpy dependency recommends building with gcc and using openblas for BLAS/LAPACK). However, when not building ``pygeosx`` other dependencies depend on python, but an existing system version works just fine, so it can be put in GEOS's ``spack.yaml`` to prevent Spack from building it. By default, Uberenv will find and use ``.uberenv_config.json`` to build GEOS, but you can use the ``--project-json`` command line option to target ``scripts/pygeosx_configs/pygeosx.json`` to build ``pygeosx``:

.. code-block:: console

    ./scripts/uberenv/uberenv.py --spack-config-dir=/path/to/your/config/directory/ --spec="%clang@14.0.6" --project-json="scripts/pygeosx_configs/pygeosx.json"

.. note::
    When building ``pygeosx``, Spack will build various python packages, however by default they are not installed in python. There are various ways of accomplishing `this <https://spack.readthedocs.io/en/latest/basic_usage.html#extensions-python-support>`_, but the recommended approach is to use spack environments. Once you build ``pygeosx`` using Uberenv, Spack will create a view that ensures the Spack-built python can find the built python packages. For example, with a default ``uberenv_libs`` directory of dependencies, the path to the view of python will be ``uberenv_libs/._view/*/bin/python3``. If you want to use your ``pygeosx`` python3 executable in GEOS, you will need to update your host-config's ``Python3_ROOT_DIR`` and ``Python3_EXECUTABLE`` to the path to Spack's view of python.

Build Configuration
-------------------

.. warning::
	The spack build system is undergoing updates. The ``petsc`` variant and others are still a work in progress.

The GEOS Spack package has a lot of options for controlling which dependencies you would like to build and how you'd like them built. The GEOS Spack package file is at ```scripts/spack_packages/packages/geosx/package.py <https://github.com/GEOS-DEV/GEOS/tree/develop/scripts/spack_packages/packages/geosx/package.py>`_.`` The variants for the package are as follows

.. literalinclude:: ../../../../scripts/spack_packages/packages/geosx/package.py
   :language: python
   :start-after: # SPHINX_BEGIN_VARIANTS
   :end-before: # SPHINX_END_VARIANTS

For example if you wanted to build with GCC 8.3.1, without Caliper and with Hypre as the Linear Algebra Interface, your spec would be ``%gcc@8.3.1 ~caliper lai=hypre``.

The GEOS Spack package lists out the libraries that GEOS depends ons. Currently these dependencies are

.. literalinclude:: ../../../../scripts/spack_packages/packages/geosx/package.py
   :language: python
   :start-after: # SPHINX_BEGIN_DEPENDS
   :end-before: # SPHINX_END_DEPENDS

Using the Spack spec syntax you can inturn specify variants for each of the dependencies of GEOS. So for example if you could modify the spec above to build RAJA in debug by using ``%gcc@8.3.1 ~caliper lai=hypre ^raja build_type=Debug``. When building with Uberenv Spack should print out a table containing the full spec for every dependency it will build. If you would like to look at the variants for say RAJA in more detail you can find the package file at ``uberenv_libs/spack/var/spack/repos/builtin/packages/raja/package.py``.

Adding a Dependency (Advanced)
------------------------------

Adding a dependency to GEOS is straight forward if the dependency already builds with Spack. If that is the case then all you need to do is add a ``depends_on('cool-new-library')`` to the GEOS ``package.py`` file. If however the dependency doesn't have a Spack package, you will have to add one by creating a ``cool-new-library/package.yaml`` file in the ``scripts/spack_packages/packages`` directory and adding the logic to build it there.

Oftentimes (unfortunately), even when a package already exists, it might not work out of the box for your system. In this case copy over the existing ``package.py`` file from the Spack repository into ``scripts/spack_packages/packages/cool-new-library/package.py``, as if you were adding a new package, and perform your modifications there. Once you have the package working, copy the package back into the Spack repository (running Uberenv should do this for you) and commit+push your changes to Spack.
