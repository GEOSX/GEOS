import os
import glob
import xml.etree.ElementTree as ElementTree
import time
import datetime
import functools
import shutil
import subprocess
import sys
import argparse


class Status:
    """ Represents the status of a job. """

    NOT_SUBMITTED = "NOT_SUBMITTED"
    SUBMITTED = "SUBMITTED"
    SUCCESS = "SUCCESS"
    FAILURE = "FAILURE"


def getMostCubeLikeRepresentation( n ):
    """
    Of all possible boxes of integer dimensions and volume n return the dimensions of the most cube-like box.

    Args:
        n: The volume of the box.

    Returns:
        The dimensions of the most cube-like box returned as a triple.
    """
    minTriple = None
    minSurfaceArea = 1e99
    for i in range( 1, n + 1 ):
        if n % i != 0:
            continue

        n_over_i = n // i
        for j in range( 1, n_over_i + 1 ):
            if n_over_i % j != 0:
                continue

            k = n_over_i // j
            surfaceArea = 2 * i * j + 2 * i * k + 2 * j * k
            if surfaceArea < minSurfaceArea:
                minSurfaceArea = surfaceArea
                minTriple = ( i, j, k )

    return minTriple


def parseListFromString( listString ):
    """
    Given a comma separated string construct a list. The list may optionally be enclosed in {}.

    Args:
        listString: The string to parse.

    Returns:
        A list of strings.
    """
    newList = listString.strip(" {}").split(",")
    return list(map( lambda x: x.strip(), newList ))


def createDirectory( dirPath, clean=False ):
    """
    Create a new directory structure, it is ok if the directory already exists.

    Args:
        dirPath: Path to the directory to create.
        clean: True if the directory is to be cleaned.
    """
    if os.path.exists( dirPath ):
        if clean:
            shutil.rmtree( dirPath )
        elif os.path.isdir( dirPath ):
            return
        else:
            raise Exception( "{} already exists but it is not a directory.".format( dirPath ) )

    os.makedirs( dirPath )


def getToolCommand( name ):
    if name.startswith("rocprof"):
        return "rocprof -i {} -o {}".format( *name.split(",")[1:] ).split(" ")
    elif name.startswith("miperf"):
        return "miperf {} -n geosx-miperf -c "
    else:
        return []

class Machine:
    """
    An class that implements an abstract computing system.

    Concrete classes are intended to derive from Machine and provide
    specific submission commands for that machine.

    Attributes:
        name: The name of the Machine.
    """

    def __init__( self, name, **kwargs ):
        """
        Initialize a Machine with the given name.

        Args:
            self: The Machine to initialize.
            name: The name of the Machine.
        """
        self.name = name
        self.submissions_command_pargs = kwargs.get("submission",[])
        self.tool_commands = getToolCommand( kwargs.get("tool","") )

    def getToolCommand( self ):
        return self.tool_commands

    def getSubmissionCommand( self, nodes, tasks, threadsPerTask, timeLimit ):
        """
        Return a list containing the command necessary to launch a job with the given configuration on this machine.

        This method needs to be overridden in the derived class.

        Args:
            self: The Machine to get the submission command for.
            nodes: The number of nodes the job will use.
            tasks: The number of tasks in the job.
            threadsPerTask: The threads per task, may be None if not specified.
            timeLimit: The time limit of the job in minutes.

        Returns:
            A list containing the submission commands.
        """
        raise NotImplementedError( "This method must be overridden in the derived class." )

    def printProgress( self ):
        """
        Print the status of the submission queue.

        This method needs to be overridden in the derived class.

        Args:
            self: The Machine to get the status of.
        """
        raise NotImplementedError( "This method must be overridden in the derived class." )

    def hasGPU( self, hipcuda ):
        """ Return True if this Machine has a CUDA GPU available."""
        return False


class SlurmMachine( Machine ):
    """ A class that implements a machine that uses Slurm. """

    def getSubmissionCommand( self, nodes, tasks, threadsPerTask, timeLimit ):
        """
        Return a list containing the command necessary to launch a job with the given configuration on this machine.

        Args:
            self: The Machine to get the submission command for.
            nodes: The number of nodes the job will use.
            tasks: The number of tasks in the job.
            threadsPerTask: The threads per task, may be None if not specified.
            timeLimit: The time limit of the job in minutes, may be None if not specified.

        Returns:
            A list containing the submission commands.
        """
        command = ["srun", "-N", nodes, "-n", tasks]
        if threadsPerTask is not None:
            command += ["-c", threadsPerTask]
        if timeLimit is not None:
            command += ["-t", timeLimit]
        command.extend( self.submissions_command_pargs )

        return command

    def hasGPU( self, hipcuda ):
        if hipcuda == "hip":
            return True
        else:
            return False

    def printProgress( self ):
        """
        Print the status of the submission queue.

        Args:
            self: The Machine to get the status of.
        """
        print( subprocess.check_output( [ "squeue", "-u", os.environ[ "USER" ] ], stderr=subprocess.STDOUT ) )


class LallocMachine( Machine ):
    """ A class that implements a machine that uses Lalloc. """

    def getSubmissionCommand( self, nodes, tasks, threadsPerTask, timeLimit ):
        """
        Return a list containing the command necessary to launch a job with the given configuration on this machine.

        Args:
            self: The Machine to get the submission command for.
            nodes: The number of nodes the job will use.
            tasks: The number of tasks in the job.
            threadsPerTask: The threads per task, may be None if not specified.
            timeLimit: The time limit of the job in minutes, may be None if not specified.

        Returns:
            A list containing the submission commands.
        """
        # TODO(corbett5): something about the bsub "-Is" argument messes up the printing while jobs are running.
        command = ["lalloc", nodes, "lrun", "-n", tasks, "-g", 1]
        if threadsPerTask is not None:
            command += ["--threadsPerTask", threadsPerTask ]
        if timeLimit is not None:
            command += ["-W", timeLimit]
        command.extend( self.submissions_command_pargs )

        return command

    def printProgress( self ):
        """
        Print the status of the submission queue.

        Args:
            self: The Machine to get the status of.
        """
        print( subprocess.check_output( [ "bjobs" ], stderr=subprocess.STDOUT ) )

    def hasGPU( self, hipcuda ):
        """ Return True."""
        if hipcuda == "cuda":
            return True
        else:
            return False


def getMachine( **kwargs ):
    """ Return a Machine object corresponding to the current machine. """
    hostName = os.environ[ "HOSTNAME" ]
    if hostName.startswith( "quartz" ):
        return SlurmMachine( "quartz", **kwargs )
    elif hostName.startswith( "login2" ):
        return SlurmMachine( "crusher", **kwargs )
    elif hostName.startswith( "lassen" ):
        return LallocMachine( "lassen", **kwargs )
    else:
        raise EnvironmentError( "Cannot determine the current machine." )


class Benchmark:
    """
    An class that implements a single benchmark run of GEOSX.

    Attributes:
        geosxPath: The path to the GEOSX executable to run.
        xmlPath: The path to the input XML file to pass to GEOSX.
        name: The name of the benchmark.
        nodes: The number of nodes to run the benchmark with.
        tasks: The number of tasks to run the benchmark with.
        threadsPerTask: The number of threads per task to run the benchmark with.
            May be None if not specified.
        timeLimit: The time limit to give the scheduler in minutes, or None if not specified.
        outputDir: The directory where the benchmark is run.
        outputFile: The path to the file containing the standard output and
            standard error from the benchmark.
        runCommand: A list of arguments which appended to the submission command provides
            the full command for running this benchmark.
        process: The subproccess associated with the benchmark.
        status: The status of the benchmark.
    """

    def __init__( self, outputDir, geosxPath, xmlPath, name, nodes, tasks, threadsPerTask, timeLimit, args, autoPartition ):
        """
        Initialize a Benchmark.

        Arguments:
            self: The Benchmark to initialize.
            outputDir: The top level directory where the benchmarks are run.
            geosxPath: The path to the GEOSX executable to run.
            xmlPath: The path to the input XML file to pass to GEOSX.
            name: The name of the benchmark.
            nodes: The number of nodes to run the benchmark with.
            tasks: The number of tasks to run the benchmark with.
            threadsPerTask: The number of threads per task to run the benchmark with.
                May be None if not specified.
            timeLimit: The time limit to give the scheduler in minutes, or None if if not specified.
            args: List of extra arguments to pass to GEOSX.
            autoPartition: If true then partition arguments are generated and passed to GEOSX.
        """
        self.geosxPath = os.path.abspath( geosxPath )
        self.xmlPath = os.path.abspath( xmlPath )
        self.name = name
        self.nodes = nodes
        self.tasks = tasks
        self.threadsPerTask = threadsPerTask
        self.timeLimit = timeLimit

        xmlName = os.path.splitext( os.path.basename( self.xmlPath ) )[ 0 ]
        self.outputDir = os.path.abspath( os.path.join( outputDir, xmlName, "{}_{}".format( self.name, self.nodes ) ) )

        self.outputFile = os.path.join( self.outputDir, "output.txt" )

        self.runCommand = [self.geosxPath, "-n", "{}/{}".format( xmlName, self.name ), "-i", self.xmlPath]

        self.runCommand += args

        if autoPartition:
            pX, pY, pZ = getMostCubeLikeRepresentation( tasks )
            self.runCommand += ["-x", pX, "-y", pY, "-z", pZ]

        self.process = None
        self.status = Status.NOT_SUBMITTED

    def submit( self, machine ):
        """
        Submit to be run on the given machine.

        Arguments:
            self: The Benchmark to submit.
            machine: The machine to run on.
        """
        createDirectory( self.outputDir, True )

        submissionCommand = machine.getSubmissionCommand( self.nodes, self.tasks, self.threadsPerTask, self.timeLimit )

        submissionCommand += machine.getToolCommand( )

        submissionCommand += self.runCommand

        if machine.hasGPU("cuda"):
            submissionCommand += ["-t", "spot,profile.mpi,profile.cuda"]
        elif machine.hasGPU("hip"):
            submissionCommand += ["-t", "spot,profile.mpi,profile.hip"]
        else:
            submissionCommand += ["-t", "spot,profile.mpi"]

        submissionCommand = list(map(str, submissionCommand ))
        with open( self.outputFile, "w" ) as outputFile:
            outputFile.write( "{}\n\n".format( " ".join( submissionCommand ) ) )
            self.process = subprocess.Popen( submissionCommand, cwd=self.outputDir, stdout=outputFile, stderr=subprocess.STDOUT )

        self.status = Status.SUBMITTED

        print( "Submitted {}".format( self ) )

    def getOutputFile( self ):
        """
        Return the path to the file containing the standard output and error of the Benchmark.

        Arguments:
            self: The Benchmark to get the output of.
        """
        return self.outputFile

    def getTimingFile( self ):
        """
        Return the path to the timing file produced by the Benchmark.

        Arguments:
            self: The Benchmark to get the timing file of.
        """
        caliFiles = glob.glob( os.path.join( self.outputDir, "*.cali" ) )
        if caliFiles:
            if len( caliFiles ) > 1:
                raise Exception( "The directory {} contains multiple .cali files!".format( self.outputDir ) )
            return caliFiles[ 0 ]
        else:
            return None

    def hasCompleted( self ):
        """
        Return True iff the Benchmark has completed.

        Arguments:
            self: The Benchmark to check.
        """
        jobExited = self.status == Status.SUBMITTED and self.process.poll() is not None
        if self.getTimingFile() is not None:
            self.status = Status.SUCCESS
            print( "Completed {}".format( self ) )
        elif jobExited:
            self.status = Status.FAILURE
            print( "Failed {}".format( self ) )

        return self.status == Status.SUCCESS or self.status == Status.FAILURE

    def kill( self ):
        """
        Kill the Benchmark.

        Arguments:
            self: The Benchmark to kill.
        """
        self.process.kill()
        self.process.wait()
        self.hasCompleted()

    def __str__( self ):
        """
        Return a string represtation of the Benchmark.

        Arguments:
            self: The Benchmark to stringify.
        """
        return "Benchmark( OutputDir: {}, nodes: {}, tasks: {}, threadsPerTask: {}, timeLimit {} )".format( self.outputDir, self.nodes, self.tasks, self.threadsPerTask, self.timeLimit )


def printProgress( machine, benchmarks, startTime ):
    """
    Print the status of the submitted benchmarks.

    Arguments:
        machine: The Machine the benchmarks were submitted on.
        benchmarks: A list of the submitted Benchmarks.
        startTime: The time of the benchmarks were submitted at.
    """
    print("")

    runTime = int( time.time() - startTime )
    print( "Running for {}".format( datetime.timedelta( seconds=runTime ) ) )

    statusCount = { Status.NOT_SUBMITTED: 0, Status.SUBMITTED: 0, Status.SUCCESS: 0, Status.FAILURE: 0 }
    for benchmark in benchmarks:
        statusCount[ benchmark.status ] += 1

    for status, count in statusCount.items():
        print( "{}: {}".format( status, count ) )

    machine.printProgress()
    sys.stdout.flush()
    print("")


def allCompleted( benchmarks ):
    """
    Return True if all the benchmarks have completed.

    Arguments:
        benchmarks: A list of the submitted Benchmarks.
    """
    return functools.reduce( lambda result, benchmark : result and benchmark.hasCompleted(), benchmarks, True )


def submitAllAndWait( machine, benchmarks, timeLimit ):
    """
    Submit all of the benchmarks and wait until the complete or the time limit expires.

    Arguments:
        machine: The Machine to submit the benchmarks to.
        benchmarks: A list of the Benchmarks to submit.
        timeLimit: The number of minutes to wait after submitting the benchmarks.
    """
    startTime = time.time()
    for benchmark in benchmarks:
        benchmark.submit( machine )

    while not allCompleted( benchmarks ):
        for _ in range( 60 ):
            time.sleep( 1 )
            if allCompleted( benchmarks ):
                break

        printProgress( machine, benchmarks, startTime )

        minutesRun = ( time.time() - startTime ) / 60
        if minutesRun > timeLimit:
            print( "The time limit of {} minutes has been exceeded, killing all unfinished jobs.".format( timeLimit ) )
            for benchmark in benchmarks:
                benchmark.kill()


def getIncludedFiles( xmlFilePath ):
    includedFiles = []
    tree = ElementTree.parse( xmlFilePath )
    includeElems = tree.findall( "./Included/File")
    for elem in includeElems:
        includedFiles.append( elem.get("name") )
    return includedFiles


def getBenchmarksFromXML( xmlFilePath, machine, outputDir, geosxPath ):
    """
    Return a list of benchmarks created for the current Machine from the given XML file.

    Arguments:
        xmlFilePath: The path to the XML file.
        machine: The Machine to run the benchmarks on.
        outputDir: The top level directory all the benchmarks are to be run in.
        geosxPath: The path to the GEOSX executable to run.
    """
    benchmarks = []
    allIncludes = getIncludedFiles( xmlFilePath )
    lastLen = 0
    while( len(allIncludes) > lastLen ):
        for include in allIncludes:
            getIncludedFiles( include )
        lastLen = len(allIncludes)
    for include in allIncludes:
        shutil.copy(include,outputDir)

    tree = ElementTree.parse( xmlFilePath )
    matchingElements = tree.findall( "./Benchmarks/{}/Run".format( machine.name ) )

    if len(matchingElements) > 0:
        mesh = tree.find("./Mesh/InternalMesh")
        if mesh is None:
            raise Exception( "Scaling runs require a pre-existing internalMesh to modify per scale!" )

    for elem in matchingElements:
        name = elem.get( "name" )
        if name is None:
            raise Exception( "Expected the XML element to have an attribute 'name'." )

        nodes = elem.get( "nodes", 1 )
        nodes = int( nodes )

        tasksPerNode = elem.get( "tasksPerNode" )
        if tasksPerNode is None:
            raise Exception( "Expected the XML element to have an attribute 'tasksPerNode'." )
        tasksPerNode = int( tasksPerNode )

        threadsPerTask = elem.get( "threadsPerTask" )
        if threadsPerTask is not None:
            threadsPerTask = int( threadsPerTask )

        timeLimit = elem.get( "timeLimit" )
        if timeLimit is not None:
            timeLimit = int( timeLimit )

        args = elem.get( "args" )
        if args is None:
            args = []
        else:
            args = args.split()

        autoPartition = elem.get( "autoPartition" )
        if autoPartition is None:
            autoPartition = False
        elif autoPartition in ("Off", "On"):
            autoPartition = autoPartition == "On"
        else:
            raise Exception( "Option for autoPartition not recognized: {}. Valid options are 'On' and 'Off'.".format( autoPartition ) )

        elementTargetCalc = { "weak" : lambda scale, target : scale * nodes * tasksPerNode * target,
                              "strong" : lambda scale, target : target }
        meshSizes = elem.get( "meshSizes" )
        scales = parseListFromString( elem.get( "scaleList", "{ 1 }" ) )
        scaling = elem.get( "scaling", "weak" )

        if meshSizes is None:
            for scale in scales:
                totalNodes = nodes * int( scale )
                totalTasks = totalNodes * tasksPerNode
                benchmarks.append( Benchmark( outputDir, geosxPath, xmlFilePath, name, totalNodes, totalTasks, threadsPerTask, timeLimit, args, autoPartition ) )
        else:
            meshSizes = parseListFromString( meshSizes )
            for meshSize in meshSizes:
                for scale in scales:
                    scale = int( scale )
                    meshSize = int( meshSize )
                    elementTarget = elementTargetCalc[ scaling ]( scale, meshSize )
                    perDim = int( elementTarget ** ( 1. / 3 ) ) + 1
                    mesh.set("nx", "{{ {} }}".format( perDim ) )
                    mesh.set("ny", "{{ {} }}".format( perDim ) )
                    mesh.set("nz", "{{ {} }}".format( perDim ) )
                    xmlFilePathTuple = os.path.splitext( os.path.basename( xmlFilePath ) )
                    newXmlFilePath = os.path.join( outputDir, xmlFilePathTuple[0] + "_{}M_elem.xml".format( elementTarget // 1000000  ) )
                    tree.write( open( newXmlFilePath, 'wb' ) )
                    totalNodes = nodes * int( scale )
                    totalTasks = totalNodes * tasksPerNode
                    benchmarks.append( Benchmark( outputDir, geosxPath, newXmlFilePath, name, totalNodes, totalTasks, threadsPerTask, timeLimit, args, autoPartition ) )

    return benchmarks


def getBenchmarksFromDirectory( benchmarkDir, machine, outputDir, geosxPath ):
    """
    Return a list of benchmarks created for the current Machine from the XML files in the given directory.

    Arguments:
        benchmarkDir: The directory containing the XML files to run.
        machine: The Machine to run the benchmarks on.
        outputDir: The top level directory all the benchmarks are to be run in.
        geosxPath: The path to the GEOSX executable to run.
    """
    benchmarks = []
    for fileName in os.listdir( benchmarkDir ):
        if fileName.endswith( ".xml" ):
            benchmarks += getBenchmarksFromXML( os.path.join( benchmarkDir, fileName ), machine, outputDir, geosxPath )

    return benchmarks


def main():
    """ Parse the command line arguments, submit all found benchmarks and wait for them to finish. """

    benchmarkDir = os.path.dirname( os.path.realpath( __file__ ) )

    timeLimit = 60
    parser = argparse.ArgumentParser()
    parser.add_argument( "geosxPath", help="The path to the GEOSX executable to benchmark.")
    parser.add_argument( "outputDirectory", help="The parent directory to run the benchmarks in.")
    parser.add_argument( "-t", "--timeLimit", type=int, help="Time limit for the entire script in minutes, the default is {}.". format( timeLimit ), default=timeLimit )
    parser.add_argument( "-o", "--timingCollectionDir", help="Directory to copy the timing files to." )
    parser.add_argument( "-e", "--errorCollectionDir", help="Directory to copy the output from any failed runs to." )
    parser.add_argument( "--tool", default = "", type=str, help="Profiler other tool to use in job invocation.")
    args, unknown = parser.parse_known_args()

    print(f"Passing unknown arguments '{unknown}' to be used by the batch submission command.")

    geosxPath = os.path.abspath( args.geosxPath )
    if not os.path.isfile( geosxPath ) or not os.access( geosxPath, os.X_OK):
        raise ValueError( "geosxPath is not an executable!" )

    outputDir = os.path.abspath( args.outputDirectory )

    timeLimit = args.timeLimit

    timingCollectionDir = args.timingCollectionDir
    if timingCollectionDir is not None:
        timingCollectionDir = os.path.abspath( timingCollectionDir )

    errorCollectionDir = args.errorCollectionDir
    if errorCollectionDir is not None:
        errorCollectionDir = os.path.abspath( errorCollectionDir )

    machine = getMachine( tool = args.tool, submission = unknown )

    benchmarks = getBenchmarksFromDirectory( benchmarkDir, machine, outputDir, geosxPath )

    print( "Benchmarking GEOSX found at {}".format( geosxPath ) )
    print( "Results will be written to {}".format( outputDir ) )
    if timingCollectionDir is not None:
        print( "Timing files will be added to {}".format( timingCollectionDir ) )
    if errorCollectionDir is not None:
        print( "Output from failed benchmarks will be put in {}".format( errorCollectionDir ) )
    print( "The time limit is {} minutes.".format( timeLimit) )
    print( "Found {} benchmarks.".format( len( benchmarks ) ) )
    print( "" )

    submitAllAndWait( machine, benchmarks, timeLimit )

    # Copy the timing files from successful benchmarks to a new directory if asked.
    if timingCollectionDir is not None:
        createDirectory( timingCollectionDir )
        for benchmark in benchmarks:
            if benchmark.status == Status.SUCCESS:
                shutil.copy2( benchmark.getTimingFile(), timingCollectionDir )

    # Copy the output from the failed benchmarks to a new directory if asked.
    if errorCollectionDir is not None:
        createDirectory( errorCollectionDir, True )
        for benchmark in benchmarks:
            if benchmark.status != Status.SUCCESS and os.path.isfile( benchmark.getOutputFile() ):
                outputFile = benchmark.getOutputFile()
                newLocation = os.path.join( errorCollectionDir, os.path.relpath( outputFile, outputDir ) )
                shutil.copy2( outputFile, newLocation )

    success = True
    for benchmark in benchmarks:
        success = success and benchmark.status == Status.SUCCESS

    return not success


if __name__ == "__main__" and not sys.flags.interactive:
    sys.exit(main())
