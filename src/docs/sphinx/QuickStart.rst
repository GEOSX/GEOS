.. _QuickStart:

###############################
Quick Start Guide
###############################

The goal of this page is to get you started as quickly as possible using GEOS.
We will walk you through downloading the source, compiling the code, and testing the installation.

Before jumping to the installation process, we want to first address some frequently asked questions we get from new users.
If you are itching to get started, feel free to jump ahead to the relevant sections.

Frequently Asked Questions
==========================

Does GEOS have a graphical user interface?:
------------------------------------------------
Given the focus on rapid development and HPC environments, GEOS does not have a graphical user interface.
This is consistent with many other high performance computing packages, but we recognize it can be a deal-breaker for certain users.
For those who can get past this failing, we promise we still have a lot to offer.
In a typical workflow, you will prepare an XML-based input file describing your problem.
You may also prepare a mesh file containing geometric and property information describing, say, a reservoir you would like to simulate.
There is no shortage of GUI tools that can help you in this model building stage.
The resulting input deck is then consumed by GEOS to run the simulation and produce results.
This may be done in a terminal of your local machine or by submitting a job to a remote server.
The resulting output files can then be visualized by any number of graphical visualization programs (typically `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ or `paraview <https://www.paraview.org/>`_).
Thus, while GEOS is GUI free, the typical workflow is not.

Do I need to be a code developer to use GEOS?:
------------------------------------------------
For the moment, most users will
need to download and compile the code from source, which we readily admit this requires
a certain level of development expertise.  We try to make this process as easy as
possible, and we are working on additional deployment options to make this process easier.
Once installed, however, our goal is to make GEOS accessible to developers and non-developers alike.
Our target audience includes engineers and scientists who want to solve tough application problems, but could care less about the insides of the tool.
For those of you who *are* interested in scientific computing, however, GEOS is an open source project and we welcome external contributions.

What are the system requirements?:
------------------------------------------------
GEOS is primarily written in C++, with a focus on standards compliance and platform-to-platform portability.
It is designed to run on everything from commodity laptops to the world's most powerful supercomputers.
We regularly test the code across a variety of operating systems and compilers.
Most of these operating systems are Linux/UNIX based (e.g. Ubuntu, CentOS, Mac OSX).
We do have developers working in Windows environments, but they use a Virtual Machine or work within a docker image rather than directly in the Windows environment.
In the instructions below, we assume you have access to fairly standard development tools.
Using advanced features of GEOS, like GPU-acceleration, will of course introduce additional hardware and software requirements.

Help, I get errors while trying to download/compile/run!:
---------------------------------------------------------

Unfortunately, no set of instructions is foolproof.
It is simply impossible to anticipate every system configuration or user.
If you run into problems during the installation, we recommend the following five-step process:

#. Take a moment to relax, and then re-read the instructions carefully.
   Perhaps you overlooked a key step?  Re-read the error message(s) closely.
   Modern compilation tools are often quite helpful in reporting exactly why things fail.

#. Type a few keywords from your error into a search engine.
   It is possible someone else out there has encountered your problem before, and a well-chosen keyword can often produce an instant solution.
   Note that when a compilation fails, you may get pages and pages of errors.  Try to identify the *first* one to occur and fix that.
   One error will often trigger subsequent errors, and looking at the *last* error on the screen may not be so helpful.

#. If you encounter problems building one of the third-party libraries we depend on, check out their support pages.
   They may be able to help you more directly than we can.

#. Still stuck? Check out our `issues tracker <https://github.com/GEOS-DEV/GEOS/issues>`_, searching current or closed issues that may address your problem.
   Perhaps someone has had an identical issue, or something close.  The issue tracker has a convenient search bar where you can search for relevant keywords.
   Remember to remove the default ``is:open`` keyword to search both open and closed issues.

#. If you have exhausted the options above, it is time to seek help from the developers.
   Post an issue on our issue tracker.
   Be specific, providing as much information as possible about your system setup and the error you are encountering.
   Please be patient in this process, as we may need to correspond a few times and ask you to run additional tests.
   Most of the time, users have a slightly unusual system configuration that we haven't encountered yet, such as an older version of a particular library.
   Other times there is a legitimate bug in GEOS to be addressed.
   Take pride in the fact that you may be saving the next user from wasted time and frustration.

Repository Organization
==============================

The source for GEOS and related tools are hosted on `Github <https://github.com>`_.
We use `Git workflows <https://git-scm.com>`_ to version control our code and manage the entire development process.
On Github, we have a `GEOS Organization <https://github.com/GEOS-DEV>`_ that hosts several related repositories.

You should sign up for a free Github account, particularly if you are interested in posting issues to our issue tracker and communicating with the developers.
The main repository of interest is obviously GEOS itself: `GEOS <https://github.com/GEOS-DEV/GEOS>`_

We also rely on two types of dependencies: first-party and third-party.
First-party dependencies are projects directly associated with the GEOS effort, but kept in separate repositories because they form stand-alone tools.
For example, there is an equation-of-state package called `PVTPackage <https://github.com/GEOS-DEV/PVTPackage>`_ or the streamlined CMake-based foundation `BLT <https://github.com/LLNL/blt>`_ .
These packages are handled as `Git Submodules <https://git-scm.com/book/en/v2/Git-Tools-Submodules>`_, which provides a transparent way of coordinating multiple code development projects.
Most users will never have to worry that these modules are in fact separate projects from GEOS.

We also rely on several open-source Third-Party Libraries (TPLs) (see `thirdPartyLibs <https://github.com/GEOS-DEV/thirdPartyLibs>`_).
These are well-respected projects developed externally to GEOS.
We have found, however, that many compilation issues stem from version incompatibilities between different packages.
To address this, we provide a mirror of these TPLs, with version combinations we know play nicely together.
We also provide a build script that conveniently and consistently builds those dependencies.

Our build system will automatically use the mirror package versions by default.
You are welcome to tune your configuration, however, to point to different versions installed on your system.
If you work on an HPC platform, for example, common packages may already be available and optimized for platform hardware.
For new users, however, it may be safer to begin with the TPL mirror.

.. note::
   If you are working on an HPC platform with several other GEOS users, we often compile the TPLs in a shared location so individual users don't have to waste their storage quota.
   Inquire with your institution's point-of-contact whether this option already exists.
   For all LLNL systems, the answer is yes.

Finally, there are also several private repositories only accessible to the core development team, which we use for behind-the-scene testing and maintenance of the code.

Username and Authentication
=============================
New users should sign up for a free `Github account <https://github.com>`_.

If you intend to develop in the GEOS codebase, you may benefit from setting up your git credentials (see :ref:`GitWorkflow`).


Download
======================

It is possible to directly download the source code as a zip file.
We strongly suggest, however, that users don't rely on this option.
Instead, most users should use Git to either *clone* or *fork* the repository.
This makes it much easier to stay up to date with the latest releases and bug fixes.
If you are not familiar with the basics of Git, `here is a helpful resource <https://git-scm.com>`_ to get you started.

The tutorial here assumes you will use a https clone with no specific credentials.
Using an ssh connection pattern requires a very slight modification.
See the **Additional Notes** at the end of this section for details.

If you do not already have Git installed on your system, you will need to install it.
We recommend using a relatively recent version of Git, as there have been some notable improvements over the past few years.
You can check if Git is already available by opening a terminal and typing

.. code-block:: sh

  git --version

You'll also need the `git-lfs <https://git-lfs.github.com/>`_ large file extension.

The first task is to clone the ``GEOS`` and ``thirdPartyLibs`` repositories.
If you do not tell it otherwise, the build system will expect the GEOS and thirdPartyLibs to be parallel to each other in the directory structure.
For example,

.. code-block:: sh

  codes/
  ├── GEOS/
  └── thirdPartyLibs/

where the toplevel ``codes`` directory can be re-named and located wherever you like.
It is possible to customize the build system to expect a different structure, but for now let us assume you take the simplest approach.

First, using a terminal, create the ``codes`` directory wherever you like.

.. code-block:: sh

  cd /insert/your/desired/path/
  mkdir codes
  cd codes

Inside this directory, we can clone the GEOS repository.
We will also use some Git commands to initialize and download the submodules (e.g. ``LvArray``).

.. code-block:: sh

   git clone https://github.com/GEOS-DEV/GEOS.git
   cd GEOS
   git lfs install
   git submodule init
   git submodule update
   cd ..

If all goes well, you should have a complete copy of the GEOS source at this point.
The most common errors people encounter here have to do with Github not recognizing their authentication settings and/or repository permissions.
See the previous section for tips on ensuring your SSH is working properly.

*Note*: Previous versions of GEOS also imported the integratedTests submodule, which is not publicly available (access is limited to the core development team).
This may cause the ``git submodule update`` command to fail.
In that case, run ``git submodule deinit integratedTests`` before ``git submodule update``.
This submodule is not required for building GEOS.

.. code-block:: sh

   cd GEOS
   git submodule update --init src/cmake/blt
   git submodule update --init src/coreComponents/LvArray
   git submodule update --init src/coreComponents/fileIO/coupling/hdf5_interface
   git submodule update --init src/coreComponents/constitutive/PVTPackage
   cd ..

Once we have grabbed GEOS, we do the same for the thirdPartyLibs repository.  From the ``codes`` directory, type

.. code-block:: sh

   git clone https://github.com/GEOS-DEV/thirdPartyLibs.git
   cd thirdPartyLibs
   git lfs install
   git pull
   git submodule init
   git submodule update
   cd ..

Again, if all goes well you should now have a copy of all necessary TPL packages.

**Additional Notes:**

#. ``git-lfs`` may not function properly (or may be very slow) if your version of git and git-lfs are not current.
If you are using an older version, you may need to add ``git lfs pull`` after ``git pull`` in the above procedures.

#. You can adapt the commands if you use an ssh connection instead.
The clone ``https://github.com/GEOS-DEV/GEOS.git`` becomes ``git clone git@github.com:GEOS-DEV/GEOS.git``.
You may also be willing to insert your credentials in the command line (less secure) ``git clone https://${USER}:${TOKEN}@github.com/GEOS-DEV/GEOS.git``.

Configuration 
=============

Before proceeding, make sure to have installed all the minimal prerequisites as described in :ref:`Prerequisites`
Note that GEOS supports a variety of parallel computing models, depending on the hardware and software environment.
Advanced users are referred to the :ref:`BuildGuide` for a discussion of the available configuration options.

Before beginning, it is a good idea to have a clear idea of the flavor and version of the build tools you are using.
If something goes wrong, the first thing the support team will ask you for is this information.

.. code-block:: sh

  cpp --version
  mpic++ --version
  cmake --version

Here, you may need to replace ``cpp`` with the full path to the C++ compiler you would like to use, depending on how your path and any aliases are configured.

Defining a Host-Config File
---------------------------

GEOS compilations are driven by a CMake ``host-config`` file, which informs the build system about the compilers you are using, where various packages reside, and what options you want to enable. 

A template for creating a simple ``host-config`` is provided in ``host-configs/quick-start-template.cmake``.

.. literalinclude:: ../../../host-configs/quick-start-template.cmake
   :language: sh

The various ``set()`` commands are used to set variables that control the build. To begin, make a copy of the template file and modify the paths according to the installation locations on your system. 

We have created a number of default host-config files for common systems. You should browse them to see if any are close to your needs:
We maintain host configuration files (ending in ``.cmake``) for HPC systems at various institutions, as well as for common personal systems. 
If you cannot find one that matches your needs, we suggest starting with one of the shorter ones and modifying it as needed. 

.. note::
   If you develop a new ``host-config`` for a particular platform that may be useful for other users, please consider sharing it with the developer team.

Compilation
===========

The configuration process for both the third-party libraries (TPLs) and GEOS is managed through a Python script called ``config-build.py``. This script simplifies and automates the setup by configuring the build and install directories and by running CMake based on the options set in the host-config file 
which is passed as a command-lne argument. The ``config-build.py`` script has several command-line options. Here, we will only use some basic options and rely on default values for many others. During this build process there wil be automatically generated build and install directories for both the TPLs and the main code,
with names consistent with the name specified in the host-config by the variable ``CONFIG_NAME``, i.e. ``build-your-platform-release`` and ``install-your-platform-release``. 

All options can be visualized by running

.. code-block:: sh

   cd thirdPartyLibs
   python scripts/config-build.py -h

.. note::

   It is strongly recommended that GEOS and TPLs be configured using the same host configuration file. Below, we assume that you keep this file in, for example, ``GEOS/host-configs/your-platform.cmake``, but the exact location is up to you.

Compiling the TPLs
-------------------

.. note::

   If you are working on an HPC system with other GEOS developers, check with them to see if the TPLs have already been compiled in a shared directory. If this is the case, you can skip ahead to just compiling the main code.
   If you are working on your own machine, you will need to configure and compile both the TPLs and the main code.

We begin by configuring the third-party libraries (TPLs) using the ``config-build.py`` script. This script sets up the build directory and runs CMake to generate the necessary build files.

.. code-block:: sh

   cd thirdPartyLibs
   python scripts/config-build.py -hc ../GEOS/host-configs/your-platform.cmake -bt Release

The TPLs will be configured in a build directory named consistently with your host configuration file, i.e., ``build-your-platform-release``.

.. code-block:: sh

   cd build-your-platform-release
   make

.. note::

   Building all of the TPLs can take quite a while, so you may want to go get a cup of coffee at this point.
   Also note that you should *not* use a parallel ``make -j N`` command to try and speed up the build time.

Compiling GEOS
-------------------

Once the TPLs have been compiler, the next step is to compile the main code. The ``config-build.py`` script is used to configure the build directory. Before running the configuration script, ensure that the path to the TPLs is correctly set in the host configuration file by setting

.. code-block:: sh

   set(GEOS_TPL_DIR "/path/to/your/TPL/installation/dir" CACHE PATH "")

If you have followed these instructions, the TPLs are installed at the default location, i.e. ``/path/to/your/TPL/thirdPartyLibs/install-your-platform-release``.

.. code-block:: sh

   cd ../../GEOS
   python scripts/config-build.py -hc host-configs/your-platform.cmake -bt Release

An alternative is to set the path ``GEOS_TPL_DIR`` via a cmake command line option, e.g.

.. code-block:: sh

   python scripts/config-build.py -hc host-configs/your-platform.cmake -bt Release -D GEOS_TPL_DIR=/full/path/to/thirdPartyLibs

.. note::

   We highly recommend using full paths, rather than relative paths, whenever possible.

Once the configuration process is completed, we proceed with the compilation of the main code and the instalation of geos.  

.. code-block:: sh

   cd build-your-platform-release
   make -j4
   make install   

The parallel ``make -j 4`` will use four processes for compilation, which can substantially speed up the build if you have a multi-processor machine.
You can adjust this value to match the number of processors available on your machine.
The ``make install`` command then installs GEOS to a default location unless otherwise specified.

If all goes well, a ``geosx`` executable should now be available

.. code-block:: sh

  GEOS/install-your-platform-release/bin/geosx

Running
=================

We can do a quick check that the geosx executable is working properly by calling the executable with our help flag

.. code-block:: sh

  ./bin/geosx --help

This should print out a brief summary of the available command line arguments:

.. code-block:: sh

    USAGE: geosx -i input.xml [options]

    Options:
    -?, --help
    -i, --input,             Input xml filename (required)
    -r, --restart,           Target restart filename
    -x, --x-partitions,      Number of partitions in the x-direction
    -y, --y-partitions,      Number of partitions in the y-direction
    -z, --z-partitions,      Number of partitions in the z-direction
    -s, --schema,            Name of the output schema
    -b, --use-nonblocking,   Use non-blocking MPI communication
    -n, --name,              Name of the problem, used for output
    -s, --suppress-pinned,   Suppress usage of pinned memory for MPI communication buffers
    -o, --output,            Directory to put the output files
    -t, --timers,            String specifying the type of timer output
    --trace-data-migration,  Trace host-device data migration
    --pause-for,             Pause geosx for a given number of seconds before starting execution

Obviously this doesn't do much interesting, but it will at least confirm that the executable runs.
In typical usage, an input XML must be provided describing the problem to be run, e.g.

.. code-block:: sh

    ./bin/geosx -i your-problem.xml

In a parallel setting, the command might look something like

.. code-block:: sh

    mpirun -np 8 ./bin/geosx -i your-problem.xml -x 2 -y 2 -z 2

Note that we provide a series of :ref:`Tutorials` to walk you through the actual usage of the code, with several input examples.
Once you are comfortable the build is working properly, we suggest new users start working through these tutorials.

Testing
=================

It is wise to run our unit test suite as an additional check that everything is working properly.
You can run them in the build folder you just created.

.. code-block:: sh

  cd GEOS/build-your-platform-release
  ctest -V

This will run a large suite of simple tests that check various components of the code.
If you have access, you may also consider running the integrated tests.
Please refer to :ref:`IntegratedTests` for further information.

.. note::
   If *all* of the unit tests fail, there is likely something wrong with your installation.
   Refer to the FAQs above for how best to proceed in this situation.
   If only a few tests fail, it is possible that your platform configuration has exposed some issue that our existing platform tests do not catch.
   If you suspect this is the case, please consider posting an issue to our issue tracker (after first checking whether other users have encountered a similar issue).
