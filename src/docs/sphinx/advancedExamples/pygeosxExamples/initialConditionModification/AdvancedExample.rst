.. _pygeosxInitialConditionModification:


#######################################################
Initial Condition Modification
#######################################################

**Objectives**

At the end of this example you will know:

  - how to modify GEOS arrays using pygeosx
  - handle parallel communication with pygeosx_tools


**Input files**

This example requires an input xml and python script located at:

.. code-block:: console

  GEOS/examples/pygeosxExamples/sedovWithStressFunction


------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

This example is derived from the sedov integrated test (`GEOS/src/coreComponents/physicsSolvers/solidMechanics/integratedTests/sedov.xml`), which looks at the propagation of elastic waves due to an initial stress field.
The pygeosx interface is used to modify the initial conditions of the problem to something of our choosing.


---------------------------------
XML Configuration
---------------------------------

As before, the basic sedov input xml file for this example requires some modification in order to work with pygeosx.
First, we use the advanced xml input features to include the base problem (this path may need to be updated, depending on where you run the problem).

.. literalinclude:: ../../../../../../examples/pygeosxExamples/sedovWithStressFunction/modified_sedov.xml
  :language: xml
  :start-after: <!-- SPHINX_SEDOV_INCLUDED -->
  :end-before: <!-- SPHINX_SEDOV_INCLUDED_END -->

Next, we add a new entry to the output block Python and an entry in the Events block.
Whenever the python event is triggered, GEOS will pause and return to the controlling python script.

.. literalinclude:: ../../../../../../examples/pygeosxExamples/sedovWithStressFunction/modified_sedov.xml
  :language: xml
  :start-after: <!-- SPHINX_SEDOV_EVENTS_OUTPUTS -->
  :end-before: <!-- SPHINX_SEDOV_EVENTS_OUTPUTS_END -->


---------------------------------
Python Script
---------------------------------

Similar to the previous example, the python script begins by importing the required packages, applying the xml preprocessor, GEOS initialization, and key search.

.. literalinclude:: ../../../../../../examples/pygeosxExamples/sedovWithStressFunction/run_sedov_problem.py
  :language: python
  :start-after: # PYGEOSX_INITIALIZATION
  :end-before: # PYGEOSX_INITIALIZATION_END

The next steps rely on a python function that we use to set stress.
The argument to this function, x, is assumed to be a numpy array of element centers:

.. literalinclude:: ../../../../../../examples/pygeosxExamples/sedovWithStressFunction/run_sedov_problem.py
  :language: python
  :start-after: # PYGEOSX_STRESS_FN
  :end-before: # PYGEOSX_STRESS_FN_END

In the following section, we zero out the initial stress and then set it based on `stress_fn`.
While doing this, we use `wrapper.print_global_value_range` to check on the process.

.. literalinclude:: ../../../../../../examples/pygeosxExamples/sedovWithStressFunction/run_sedov_problem.py
  :language: python
  :start-after: # PYGEOSX_STRESS_IC
  :end-before: # PYGEOSX_STRESS_IC_END

Finally, we run the simulation.
As an optional step, we extract numpy arrays from the datastructure using different parallel approaches:

.. literalinclude:: ../../../../../../examples/pygeosxExamples/sedovWithStressFunction/run_sedov_problem.py
  :language: python
  :start-after: # PYGEOSX_MAIN_LOOP
  :end-before: # PYGEOSX_MAIN_LOOP_END


---------------------------------
Running the Problem
---------------------------------

To run the problem, you must use the specific version of python where pygeosx is installed.
This is likeley located here:

.. code-block:: console

   GEOS/[build_dir]/lib/PYGEOSX/bin/python

Note that you may need to manually install the pygeosx_tools package (and its pre-requisites) into this python distribution.
To do so, you can run the following:

.. code-block:: console

   cd GEOS/[build_dir]/lib/PYGEOSX/bin
   pip install --upgrade ../../../../src/coreComponents/python/modules/pygeosx_tools_package/

To run the code, you will call the pygeosx run script with python, and supply the typical geosx command-line arguments and any parallel arguments.
For example:

.. code-block:: console

   # Load the correct python environment
   # If you are not using a bash shell, you may need to target one of
   # the other activation scripts
   source GEOS/[build_dir]/lib/PYGEOSX/bin/activate

   # Move to the correct directory and run
   cd /path/to/problem
   python run_sedov_problem.py -i modified_sedov.xml -o results


------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.

**For more details**

  - More on advanced xml features, please see :ref:`AdvancedXMLFeatures`.
  - More on the pygeosx interface, please see :ref:`pygeosxInterface`.
