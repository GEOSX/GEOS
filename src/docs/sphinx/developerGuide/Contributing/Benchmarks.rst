.. _benchmarks:

Benchmarks
##########

In addition to the integrated tests which track code correctness we have a suite of benchmarks that track performance.


Running the benchmarks
----------------------

Because performance is system specific we currently only support running the benchmarks on the LLNL machines Quartz and Lassen. If you are on either of these machines the script ``benchmarks/runBenchmarks.py`` can be used to run the benchmarks.

::

    > python ../benchmarks/runBenchmarks.py --help
    usage: runBenchmarks.py [-h] [-t TIMELIMIT] [-o TIMINGCOLLECTIONDIR]
                            [-e ERRORCOLLECTIONDIR]
                            geosxPath outputDirectory

    positional arguments:
      geosxPath             The path to the GEOS executable to benchmark.
      outputDirectory       The parent directory to run the benchmarks in.

    optional arguments:
      -h, --help            show this help message and exit
      -t TIMELIMIT, --timeLimit TIMELIMIT
                            Time limit for the entire script in minutes, the
                            default is 60.
      -o TIMINGCOLLECTIONDIR, --timingCollectionDir TIMINGCOLLECTIONDIR
                            Directory to copy the timing files to.
      -e ERRORCOLLECTIONDIR, --errorCollectionDir ERRORCOLLECTIONDIR
                            Directory to copy the output from any failed runs to.

At a minimum you need to pass the script the path to the GEOS executable and a directory to run the benchmarks in. This directory will be created if it doesn't exist. The script will collect a list of benchmarks to be run and submit a job to the system's scheduler for each benchmark. This means that you don't need to be in an allocation to run the benchmarks. Note that this is different from the integrated tests where you need to already be in an allocation and an internal scheduler is used to run the individual tests. Since a benchmark is a measure of performance to get consistent results it is important that each time a benchmark is run it has access to the same resources. Using the system scheduler guarantees this.

In addition to whatever outputs the input would normally produce (plot files, restart files, ...) each benchmark will produce an output file ``output.txt`` containing the standard output and standard error of the run and a ``.cali`` file containing the Caliper timing data in a format that Spot_ can read.

.. note::
  A future version of the script will be able to run only a subset of the benchmarks.

Specifying a benchmark
----------------------

A group of benchmarks is specified with a standard GEOS input XML file with an extra ``Benchmarks`` block added at the top level. This block is ignored by GEOS itself and only used by the ``runBenchmarks.py`` script.

.. literalinclude:: ../../../../../inputFiles/solidMechanics/benchmarks/SSLE-io.xml
   :language: xml
   :start-after: <Problem>
   :end-before: <Solvers>

*[Source: benchmarks/SSLE-small.xml]*

The ``Benchmarks`` block consists of a block for each machine the benchmarks are to run on. Currently the only options are ``quartz``, ``lassen``, and ``crusher``.


The ``Run`` block
^^^^^^^^^^^^^^^^^

Each machine block contains a number of ``Run`` blocks each of which specify a family of benchmarks to run. Each ``Run`` block must have the following required attributes

  - ``name``: The name of the family of benchmarks, must be unique among all the other ``Run`` blocks on that system.
  - ``nodes``: An integer which specifies the base number of nodes to run the benchmark with.
  - ``tasksPerNode``: An integer that specifies the number of tasks to launch per node.

Each ``Run`` block may contain the following optional attributes

  - ``threadsPerTask``: An integer specifying the number of threads to allocate each task.
  - ``timeLimit``: An integer specifying the time limit in minutes to pass to the system scheduler when submitting the benchmark.
  - ``args``: containing any extra command line arguments to pass to GEOS.
  - ``autoPartition``: Either ``On`` or ``Off``, not specifying ``autoPartition`` is equivalent to ``autoPartition="Off"``. When auto partitioning is enabled the script will compute the number of ``x``, ``y`` and ``z`` partitions such that the the resulting partition is close to a perfect cube as possible, ie with 27 tasks ``x = 3, y = 3, z = 3`` and with 36 tasks ``x = 4, y = 3, z = 3``. This is optimal when the domain itself is a cube, but will be suboptimal otherwise.
  - ``strongScaling``: A list of unique integers specifying the factors to scale the number of nodes by. If ``N`` number are provided then ``N`` benchmarks are run and benchmark ``i`` uses ``nodes * strongScaling[ i ]`` nodes. Not specifying ``strongScaling`` is equivalent to ``strongScaling="{ 1 }"``.

Looking at the example ``Benchmarks`` block above on Lassen one benchmark from the ``OMP_CUDA`` family will be run with one node and one task. Four benchmarks from the ``MPI_OMP_CUDA`` family will be run with one, two, four and eight nodes and four tasks per node.

Note that specifying a time limit for each benchmark family can greatly decrease the time spent waiting in the scheduler's queue. A good rule of thumb is that the time limit should be twice as long as it takes to run the longest benchmark in the family.


Adding a benchmark problem
---------------------------

To add a new group of benchmarks you need to create an XML input file describing the problem to be run. Then you need to add the ``Benchmarks`` block described above which specifies the specific benchmarks. Finally add a symbolic link to the input file in ``benchmarks`` and run the benchmarks to make sure everything works as expected.


Viewing the results
-------------------

Each night the NightlyTests_ repository runs the benchmarks on both Quartz and Lassen, the ``timingFiles`` directory contains all of the resulting caliper output files. If you're on LC then these files are duplicated at ``/usr/gapps/GEOSX/timingFiles/`` and if you have LC access you can view them in Spot_. You can also open these files in Python and analyse them (See :ref:`opening-spot-caliper-files-in-python`).

If you want to run the benchmarks on your local branch and compare the results with develop you can use the ``benchmarks/compareBenchmarks.py`` python script. This requires that you run the benchmarks on your branch and on develop. It will print out a table with the initialization time speed up and run time speed up, so a run speed up of of 2x means your branch runs twice as fast as develop where as a initialization speed up of 0.5x means the set up takes twice as long.

.. note::
  A future version of the script will be able to pull timing results straight from the ``.cali`` files so that if you have access to the NightlyTests_ timing files you won't need to run the benchmarks on develop. Furthermore it will be able to provide more detailed information than just initialization and run times.

.. _NightlyTests: https://github.com/GEOS-DEV/NightlyTests
.. _Spot: https://lc.llnl.gov/spot2/?sf=/usr/gapps/GEOSX/timingFiles
