.. _pygeosxInSituDataMonitor:


#######################################################
In Situ Data Monitor
#######################################################

**Objectives**

At the end of this example you will know:

  - how to run a problem using the pygeosx interface,
  - how to process advanced xml features using pygeosx,
  - how to extract and monitor values within the GEOS datastructure in real-time


**Input files**

This example requires two input xml files and one python script located at:

.. code-block:: console

  GEOS/examples/pygeosxExamples/hydraulicFractureWithMonitor


------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

This example is derived from this basic example: :ref:`TutorialHydraulicFractureWithAdvancedXML`, which solves for the propagation of a single hydraulic fracture within a heterogeneous reservoir.
The pygeosx interface is used to monitor the maximum hydraulic aperture and fracture extents over time.


---------------------------------
XML Configuration
---------------------------------

The input xml file for this example requires some modification in order to work with pygeosx.
First, we use the advanced xml input features to include the base problem and override the `table_root` parameter that points to the table files.
Note that these paths will need to be updated if you run this problem outside of the example directory.

.. literalinclude:: ../../../../../../examples/pygeosxExamples/hydraulicFractureWithMonitor/hydraulicFracture.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACTURE_INCLUDED_PARAMETERS -->
  :end-before: <!-- SPHINX_HYDROFRACTURE_INCLUDED_PARAMETERS_END -->

Next, we add a new entry to the output block Python and an entry in the Events block.
Whenever the python event is triggered, GEOS will pause and return to the controlling python script (in this case, every 10 cycles).

.. literalinclude:: ../../../../../../examples/pygeosxExamples/hydraulicFractureWithMonitor/hydraulicFracture.xml
  :language: xml
  :start-after: <!-- SPHINX_HYDROFRACTURE_EVENTS_OUTPUTS -->
  :end-before: <!-- SPHINX_HYDROFRACTURE_EVENTS_OUTPUTS_END -->


---------------------------------
Python Script
---------------------------------

Problems that use the pygeosx interface are driven by a custom python script.
To begin, we import a number of packages and check whether this is a parallel run.
The custom packages include pygeosx, which provides an interface to GEOS, and pygeosx_tools, which provides a number of common tools for working with the datastructure and dealing with parallel communication.

.. literalinclude:: ../../../../../../examples/pygeosxExamples/hydraulicFractureWithMonitor/hydraulicFractureWithMonitor.py
  :language: python
  :start-after: # PYGEOSX_SETUP
  :end-before: # PYGEOSX_SETUP_END

In the next step, we apply the xml preprocessor to resolve the advanced xml features.
Note that this step will modify the input arguments to reflect the location of the compiled xml file, which is processed directly by GEOS.
The script then initializes GEOS and receives the `problem` handle, which is the scripts view into the datastructure.
There is an opportunity to interact with the GEOS before the initial conditions are set, which we do not use in this example.

.. literalinclude:: ../../../../../../examples/pygeosxExamples/hydraulicFractureWithMonitor/hydraulicFractureWithMonitor.py
  :language: python
  :start-after: # PYGEOSX_INITIALIZATION
  :end-before: # PYGEOSX_INITIALIZATION_END

To extract information from the problem, you need to know the full path (or 'key') to the target object.
These keys can be quite long, and can change depending on the xml input.
In the next step, we use a method from the pygeosx_tools package to search for these keys using a list of keywords.
If the keys are known beforehand, then this step could be skipped.
Note that these functions will throw an error if they do not find a matching key, or if they find multiple matching keys.

.. literalinclude:: ../../../../../../examples/pygeosxExamples/hydraulicFractureWithMonitor/hydraulicFractureWithMonitor.py
  :language: python
  :start-after: # PYGEOSX_KEY_SEARCH
  :end-before: # PYGEOSX_KEY_SEARCH_END

Next, we setup a dictionary that will allow us to use pygeosx_tools to automatically query the problem.
The root level of this dictionary contains the target keys (fracture location and aperture) and the required `time` key.
These each point to a sub-dictionary that holds an axis label, a scale factor, and an empty list to hold the time history.
The target dictionaries also hold an entry `fhandle`, which contains a matplotlib figure handle that we can use to display the results.

.. literalinclude:: ../../../../../../examples/pygeosxExamples/hydraulicFractureWithMonitor/hydraulicFractureWithMonitor.py
  :language: python
  :start-after: # PYGEOSX_QUERY_SETUP
  :end-before: # PYGEOSX_QUERY_SETUP_END

After setting up the problem, we enter the main problem loop.
Upon calling `pygeosx.run()`, the code will execute until a Python event is triggered in the Event loop.
At those points, we have the option to interact with the problem before continuing processing.
Here, we use pygeosx_tools to query the datastructure and occasionaly plot the results to the screen.

.. literalinclude:: ../../../../../../examples/pygeosxExamples/hydraulicFractureWithMonitor/hydraulicFractureWithMonitor.py
  :language: python
  :start-after: # PYGEOSX_QUERY_SETUP
  :end-before: # PYGEOSX_QUERY_SETUP_END


---------------------------------
Manual Query
---------------------------------

To obtain and manually inspect an object in the problem, you can use the methods in `pygeosx_tools.wrapper`.
These are designed to handle any parallel communication that may be required in your analysis.
For example, to get the fracture aperture as a numpy array, you could call:

.. code-block:: python

  from pygeosx_tools import wrapper

  # (problem initialization / configuration)

  # Grab aperture as a numpy array, using three different approaches
  
  # Local copy (the write flag indicates that we do not plan to modify the result)
  aperture_local = wrapper.get_wrapper(problem, fracture_aperture_key, write_flag=False)

  # Global copy on the root rank
  aperture_global = wrapper.gather_wrapper(problem, fracture_aperture_key)

  # Global copy on the all ranks
  aperture_global = wrapper.allgather_wrapper(problem, fracture_aperture_key)


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
   srun -n 36 -ppdebug python hydraulicFractureWithMonitor.py -i hydraulicFracture.xml -x 6 -y 2 -z 3 -o hf_results


------------------------------------------------------------------
To go further
------------------------------------------------------------------

**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.

**For more details**

  - More on advanced xml features, please see :ref:`AdvancedXMLFeatures`.
  - More on the pygeosx interface, please see :ref:`pygeosxInterface`.
