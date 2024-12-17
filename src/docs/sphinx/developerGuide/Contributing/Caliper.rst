*****************************
Basic profiling with CALIPER
*****************************

GEOS is equipped with `Caliper <https://github.com/LLNL/Caliper>`_ timers.
We integrate Caliper into GEOS by marking source-code sections of interest such as computational kernels or initialization steps.
Caliper is included in the GEOS TPL library and is built by adding the following cmake configuration to a host-config file.

.. code-block:: sh

   option( ENABLE_CALIPER "Enables CALIPER" On )


GEOS/Caliper Annotation Macros
=====================================

The following macros may be used to annotate GEOS:

* ``GEOS_MARK_SCOPE(name)`` - Marks a scope with the given name.
* ``GEOS_MARK_FUNCTION`` - Marks a function with the name of the function. The name includes the namespace the function is in but not any of the template arguments or parameters. Therefore overloaded function all show up as one entry. If you would like to mark up a specific overload use ``GEOS_MARK_SCOPE`` with a unique name. 
* ``GEOS_MARK_BEGIN(name)`` - Marks the beginning of a user defined code region. 
* ``GEOS_MARK_END(name)`` - Marks the end of user defined code region.

The 2 first macros also generate annotations for NVTX is ENABLE_CUDA_NVTOOLSEXT is activated through CMake.

Configuring Caliper
=================================
  
Caliper configuration is done by specifying a string to initialize Caliper with via the
`-t` option. A few options are listed below but we refer the reader to
`Caliper Config <https://software.llnl.gov/Caliper/BuiltinConfigurations.html>`_ for the full Caliper tutorial.

* ``-t runtime-report,max_column_width=200`` Will make Caliper print aggregated timing information to standard out, with a column width large enought that it doesn't truncate most function names.
* ``-t runtime-report,max_column_width=200,profile.cuda`` Does the same as the above, but also instruments CUDA API calls. This is only an option when building with CUDA.
* ``-t runtime-report,aggregate_across_ranks=false`` Will make Caliper write per rank timing information to standard out. This isn't useful when using more than one rank but it does provide more information for single rank runs.
* ``-t spot()`` Will make Caliper output a `.cali` timing file that can be viewed in the Spot web server.


Using Adiak
=================================
Adiak is a library that allows the addition of meta-data to the Caliper Spot output, it is enabled with Caliper.
This meta-data allows you to easily slice and dice the timings available in the Spot web server. To export meta-data
use the `adiak::value` function.

See `Adiak API <https://github.com/LLNL/Adiak/blob/f27ba674b88c2435e5e3245acbda9fc0a57bf88f/docs/Adiak%20API.docx>`_
for the full Adiak documentation.


Using Spot
=================================
To use Spot you will need an LC account and a directory full of `.cali` files you would like to analyse.
Point your browser to `Spot <https://lc.llnl.gov/spot2>`_ and open up the directory containing the timing files.

.. _opening-spot-caliper-files-in-python:

Opening Spot caliper files in Python
====================================

An example Python program for analyzing Spot Caliper files in Python is provided below. Note that it requires ``pandas`` and ``hatchet`` both of which can be installed with a package manager. In addition it requires that ``cali-query`` is in the ``PATH`` variable, this is built with Caliper so we can just point it into the TPLs.

.. code-block:: Python

  import sys
  import subprocess
  import json
  import os

  import pandas as pd
  from IPython.display import display, HTML

  # Import hatchet, on LC this can be done by adding hatchet to PYTHONPATH
  sys.path.append('/usr/gapps/spot/live/hatchet')
  import hatchet as ht

  # Add cali-query to PATH
  cali_query_path = "/usr/gapps/GEOSX/thirdPartyLibs/2020-06-12/install-quartz-gcc@8.1.0-release/caliper/bin"
  os.environ["PATH"] += os.pathsep + cali_query_path

  CALI_FILES = [ 
  { "cali_file": "/usr/gapps/GEOSX/timingFiles/200612-04342891243.cali", "metric_name": "avg#inclusive#sum#time.duration"}, 
  { "cali_file": "/usr/gapps/GEOSX/timingFiles/200611-044740108300.cali", "metric_name": "avg#inclusive#sum#time.duration"}, 
  ]

  grouping_attribute = "prop:nested"
  default_metric = "avg#inclusive#sum#time.duration" 
  query = "select %s,sum(%s) group by %s format json-split" % (grouping_attribute, default_metric, grouping_attribute)

  gf1 = ht.GraphFrame.from_caliper(CALI_FILES[0]['cali_file'], query)
  gf2 = ht.GraphFrame.from_caliper(CALI_FILES[1]['cali_file'], query)

  # Print the tree representation using the default metric
  # Also print the resulting dataframe with metadata
  print(gf1.tree(color=True, metric="sum#"+default_metric))
  display(HTML(gf1.dataframe.to_html()))

  # Print the tree representation using the default metric
  # Also print the resulting dataframe with metadata
  print(gf2.tree(color=True, metric="sum#"+default_metric))
  display(HTML(gf2.dataframe.to_html()))

  # Compute the speedup between the first two cali files (exlusive and inclusive metrics only)
  gf3 = (gf1 - gf2) / gf2
  print(gf3.tree(color=True, metric="sum#"+default_metric))

  # Compute the difference between the first two cali files (exclusive and inclusive metrics only)
  # Print the resulting tree
  gf4 = gf1 - gf2
  print(gf4.tree(color=True, metric="sum#"+default_metric))

  # Compute the sum of the first two cali files (exclusive and inclusive metrics only)
  # Print the resulting tree
  gf5 = gf1 + gf2
  print(gf5.tree(color=True, metric="sum#"+default_metric))
