PVT Driver
===============

.. contents:: Table of Contents
    :depth: 3

Introduction
------------

When calibrating fluid material parameters to experimental or other reference data, it can be a hassle to launch a full flow simulation just to confirm density, viscosity, and other fluid properties are behaving as expected.  
Instead, GEOS provides a ``PVTDriver`` allowing the user to test fluid property models for a well defined set of pressure, temperature, and composition conditions.  
The driver itself is launched like any other GEOS simulation, but with a particular XML structure:

.. code-block:: sh

   ./bin/geosx -i myFluidTest.xml

This driver will work for any multi-phase fluid model (e.g. black-oil, co2-brine, compositional multiphase) enabled within GEOS.    

XML Structure
-------------
A typical XML file to run the driver will have several key elements. 
Here, we will walk through an example file included in the source tree at

.. code-block:: sh

   src/coreComponents/unitTests/constitutiveTests/testPVT_docExample.xml

The first thing to note is that the XML file structure is identical to a standard GEOS input deck.  
In fact, once the constitutive block is calibrated, one could start adding solver and discretization blocks to the same file to create a proper field simulation.  
This makes it easy to go back and forth between calibration and simulation.

The first step is to define a parameterized fluid model to test.
Here, we create a particular type of CO2-Brine mixture:

.. literalinclude:: ../../unitTests/constitutiveTests/testPVT_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_PVTDRIVER_CONSTITUTIVE_START --> 
  :end-before: <!-- SPHINX_PVTDRIVER_CONSTITUTIVE_END -->

We also define two time-history functions for the pressure (Pascal units) and temperature (Kelvin units) conditions we want to explore.

.. literalinclude:: ../../unitTests/constitutiveTests/testPVT_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_PVTDRIVER_FUNCTIONS_START --> 
  :end-before: <!-- SPHINX_PVTDRIVER_FUNCTIONS_END -->

Note that the time-axis here is just a pseudo-time, allowing us to `parameterize <https://en.wikipedia.org/wiki/Parametric_equation>`_ arbitrarily complicated paths through a (pressure,temperature) diagram.
The actual time values have no impact on the resulting fluid properties.
Here, we fix the temperature at 350K and simply ramp up pressure from 1 MPa to 50 MPa:

A ``PVTDriver`` is then added as a ``Task``, a particular type of executable event often used for simple actions. 

.. literalinclude:: ../../unitTests/constitutiveTests/testPVT_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_PVTDRIVER_TASKS_START --> 
  :end-before: <!-- SPHINX_PVTDRIVER_TASKS_END -->

The driver itself takes as input the fluid model, the pressure and temperature control functions, and a "feed composition." 
The latter is the mole fraction of each component in the mixture to be tested.
The ``steps`` parameter controls how many steps are taken along the parametric (P,T) path.
Results will be written in a simple ASCII table format (described below) to the file ``output``.
The ``logLevel`` parameter controls the verbosity of log output during execution. 
 
The driver task is added as a ``SoloEvent`` to the event queue.  
This leads to a trivial event queue, since all we do is launch the driver and then quit.

.. literalinclude:: ../../unitTests/constitutiveTests/testPVT_docExample.xml
  :language: xml
  :start-after: <!-- SPHINX_PVTDRIVER_EVENTS_START --> 
  :end-before: <!-- SPHINX_PVTDRIVER_EVENTS_END -->

Internally, the driver uses a simple form of time-stepping to advance through the (P,T) steps. 
This timestepping is handled independently of the more complicated time-stepping pattern used by physics ``Solvers`` and coordinated by the ``EventManager``.  
In particular, in the XML file above, the ``maxTime`` parameter in the ``Events`` block is an event manager control, controlling when/if certain events occur.  
Once launched, the PVTDriver internally determines its own max time and timestep size using a combination of the input functions' time coordinates and the requested number of loadsteps.  
It is therefore helpful to think of the driver as an instantaneous *event* (from the event manager's point of view), but one which has a separate, internal clock.

Parameters
----------
The key XML parameters for the PVTDriver are summarized in the following table:

.. include:: /coreComponents/schema/docs/PVTDriver.rst

Output Format
-------------
The ``output`` key is used to identify a file to which the results of the simulation are written.  
If this key is omitted, or the user specifies ``output="none"``, file output will be suppressed.  
The file is a simple ASCII format with a brief header followed by test data:

.. code:: sh

  # column 1 = time
  # column 2 = pressure
  # column 3 = temperature
  # column 4 = density
  # columns 5-6 = phase fractions
  # columns 7-8 = phase densities
  # columns 9-10 = phase viscosities
  0.0000e+00 1.0000e+06 3.5000e+02 1.5581e+01 1.0000e+00 4.1138e-11 1.5581e+01 1.0033e+03 1.7476e-05 9.9525e-04
  2.0408e-02 2.0000e+06 3.5000e+02 3.2165e+01 1.0000e+00 4.1359e-11 3.2165e+01 1.0050e+03 1.7601e-05 9.9525e-04
  4.0816e-02 3.0000e+06 3.5000e+02 4.9901e+01 1.0000e+00 4.1563e-11 4.9901e+01 1.0066e+03 1.7778e-05 9.9525e-04
  ...

Note that the number of columns will depend on how many phases and components are present.
In this case, we have a two-phase, two-component mixture.
The total density is reported in column 4, while phase fractions, phase densities, and phase viscosities are reported in subsequent columns.
If the ``outputCompressibility`` flag is activated, an extra column will be added for the total fluid compressibility after the density.
This is defined as :math:`c_t=\frac{1}{\rho_t}\left(\partial{\rho_t}/\partial P\right)` where :math:`\rho_t` is the total density.
If the ``outputMassDensity`` flag is activated, extra columns will be added for the mass density of each phase.
The number of columns will also depend on whether the ``outputPhaseComposition`` flag is activated or not. If it is activated, there will be an extra column for the mole fraction of each component in each phase.
The phase order will match the one defined in the input XML (here, the co2-rich phase followed by the water-rich phase).
This file can be readily plotted using any number of plotting tools.  Each row corresponds to one timestep of the driver, starting from initial conditions in the first row.

Unit Testing
------------

The development team also uses the PVTDriver to perform unit testing on the various fluid models within GEOS.  
The optional argument ``baseline`` can be used to point to a previous output file that has been validated  (e.g. against experimental benchmarks or reference codes).  
If such a file is specified, the driver will perform a testing run and then compare the new results against the baseline.  
In this way, any regressions in the fluid models can be quickly identified.

Developers of new models are encouraged to add their own baselines to ``src/coreComponents/constitutive/unitTests``. 
Adding additional tests is straightforward:

1. Create a new xml file for your test in ``src/coreComponents/constitutive/unitTests`` or (easier) add extra blocks to the existing XML at ``src/coreComponents/constitutive/unitTests/testPVT.xml``.  
For new XMLs, we suggest using the naming convention ``testPVT_myTest.xml``, so that all tests will be grouped together alphabetically.  
Set the ``output`` file to ``testPVT_myTest.txt``, and run your test.  
Validate the results however is appropriate.
If you have reference data available for this validation, we suggest archiving it in the ``testPVT_data/`` subdirectory, with a description of the source and formatting in the file header.
Several reference datasets are included already as examples.
This directory is also a convenient place to store auxiliary input files like PVT tables.

2. This output file will now become your new baseline.  
Replace the ``output`` key with ``baseline`` so that the driver can read in your file as a baseline for comparison.  
Make sure there is no remaining ``output`` key, or set ``output=none``, to suppress further file output.  
While you can certainly write a new output for debugging purposes, during our automated unit tests we prefer to suppress file output.  
Re-run the driver to confirm that the comparison test passes.

3. Modify ``src/coreComponents/constitutive/unitTests/CMakeLists.txt`` to enable your new test in the unit test suite.  
In particular, you will need to add your new XML file to the existing list in the ``gtest_pvt_xmls`` variable.  
Note that if you simply added your test to the existing ``testPVT.xml`` file, no changes are needed.

.. code:: sh

  set( gtest_pvt_xmls
       testPVT.xml
       testPVT_myTest.xml
     )

4. Run ``make`` in your build directory to make sure the CMake syntax is correct

5. Run ``ctest -V -R PVT`` to run the PVT unit tests.  Confirm your test is included and passes properly.

If you run into troubles, do not hesitate to contact the development team for help.
