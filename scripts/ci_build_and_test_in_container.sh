#!/bin/bash
set -o pipefail

export PYTHONDONTWRITEBYTECODE=1

printenv

SCRIPT_NAME=$0
echo "Running CLI ${SCRIPT_NAME} $@"


# docs.docker.com/config/containers/resource_constraints
# Inside the container, tools like free report the host's available swap, not what's available inside the container.
# Don't rely on the output of free or similar tools to determine whether swap is present.
echo "running free -g"
free -g

# The or_die function run the passed command line and
# exits the program in case of non zero error code
function or_die () {
    "$@"
    local status=$?

    if [[ $status != 0 ]] ; then
        echo ERROR $status command: $@
        exit $status
    fi
}

function usage () {
>&2 cat << EOF
Usage: $0
  --build-exe-only
      Request for the build of geos only.
  --cmake-build-type ...
      One of Debug, Release, RelWithDebInfo and MinSizeRel. Forwarded to CMAKE_BUILD_TYPE.
  --code-coverage
      run a code build and test.
  --data-basename output.tar.gz
      If some data needs to be extracted from the build, the argument will define the tarball. Has to be a `tar.gz`.
  --exchange-dir /path/to/exchange
      Folder to share data with outside of the container.
  --host-config host-config/my_config.cmake
      The host-config. Path is relative to the root of the repository.
  --install-dir-basename GEOS-e42ffc1
      GEOS installation basename.
  --no-install-schema
      Do not install the xsd schema.
  --no-run-unit-tests
      Do not run the unit tests (but they will be built).
  --nproc N
      Number of cores to use for the build.
  --repository /path/to/repository
      Internal mountpoint where the geos repository will be available. 
  --run-integrated-tests
      Run the integrated tests. Then bundle and send the results to the cloud.
  --sccache-credentials credentials.json
      Basename of the json credentials file to connect to the sccache cloud cache.
  --test-code-style
  --test-documentation
  -h | --help
EOF
exit 1
}

# First working in the root of the cloned repository.
# Then we'll move to the build dir.
or_die cd $(dirname $0)/..

# Parsing using getopt
args=$(or_die getopt -a -o h --long build-exe-only,cmake-build-type:,code-coverage,data-basename:,exchange-dir:,host-config:,install-dir-basename:,no-install-schema,no-run-unit-tests,nproc:,repository:,run-integrated-tests,sccache-credentials:,test-code-style,test-documentation,help -- "$@")

# Variables with default values
BUILD_EXE_ONLY=false
GEOSX_INSTALL_SCHEMA=true
HOST_CONFIG="host-configs/environment.cmake"
RUN_UNIT_TESTS=true
RUN_INTEGRATED_TESTS=false
UPLOAD_TEST_BASELINES=false
TEST_CODE_STYLE=false
TEST_DOCUMENTATION=false
CODE_COVERAGE=false
NPROC="$(nproc)"

eval set -- ${args}
while :
do
  case $1 in
    --build-exe-only)
      BUILD_EXE_ONLY=true
      RUN_UNIT_TESTS=false
      shift;;
    --cmake-build-type)      CMAKE_BUILD_TYPE=$2;        shift 2;;
    --data-basename)
      DATA_BASENAME=$2
      DATA_BASENAME_WE=${DATA_BASENAME%%.*}
      DATA_BASENAME_EXT=${DATA_BASENAME#*.}
      if [[ ${DATA_BASENAME_EXT} != "tar.gz" ]] ; then
          echo "The script ${SCRIPT_NAME} can only pack data into a '.tar.gz' file."
          exit 1
      fi
      unset DATA_BASENAME DATA_BASENAME_EXT
      shift 2;;
    --exchange-dir)          DATA_EXCHANGE_DIR=$2;       shift 2;;
    --host-config)           HOST_CONFIG=$2;             shift 2;;
    --install-dir-basename)  GEOSX_DIR=${GEOSX_TPL_DIR}/../$2; shift 2;;
    --no-install-schema)     GEOSX_INSTALL_SCHEMA=false; shift;;
    --no-run-unit-tests)     RUN_UNIT_TESTS=false;       shift;;
    --nproc)                 NPROC=$2;                   shift 2;;
    --repository)            GEOS_SRC_DIR=$2;            shift 2;;
    --run-integrated-tests)  RUN_INTEGRATED_TESTS=true;  shift;;
    --upload-test-baselines) UPLOAD_TEST_BASELINES=true; shift;;
    --code-coverage)         CODE_COVERAGE=true;         shift;;
    --sccache-credentials)   SCCACHE_CREDS=$2;           shift 2;;
    --test-code-style)       TEST_CODE_STYLE=true;       shift;;
    --test-documentation)    TEST_DOCUMENTATION=true;    shift;;
    -h | --help)             usage;                      shift;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break;;
    *) >&2 echo Unsupported option: $1
       usage;;
  esac
done

if [[ -z "${GEOS_SRC_DIR}" ]]; then
  echo "Variable GEOS_SRC_DIR is either empty or not defined. Please define it using '--repository'."
  exit 1
fi

if [[ -z "${GEOSX_DIR}" ]]; then
  echo "Installation folder undefined. Set to default value '/dev/null'. You can define it using '--install-dir-basename'."
  GEOSX_DIR=/dev/null
fi

if [[ ! -z "${SCCACHE_CREDS}" ]]; then
  # The credential json file is available at the root of the geos repository.
  # We hereafter create the config file that points to it.
  # We use this file since it's managed by the 'google-github-actions/auth' actions.
  or_die mkdir -p ${HOME}/.config/sccache
  or_die cat <<EOT >> ${HOME}/.config/sccache/config
[cache.gcs]
rw_mode = "READ_WRITE"
cred_path = "${GEOS_SRC_DIR}/${SCCACHE_CREDS}"
bucket = "geos-dev"
key_prefix = "sccache"
EOT

  # To use `sccache`, it's enough to tell `cmake` to launch the compilation using `sccache`.
  # The path to the `sccache` executable is available through the SCCACHE environment variable.
  SCCACHE_CMAKE_ARGS="-DCMAKE_CXX_COMPILER_LAUNCHER=${SCCACHE} -DCMAKE_CUDA_COMPILER_LAUNCHER=${SCCACHE}"

  if [ -n "${DOCKER_CERTS_DIR}" ] && [ -n "${DOCKER_CERTS_UPDATE_COMMAND}" ]; then
    echo "updating certificates."
    for file in "${DOCKER_CERTS_DIR}"/llnl/*.crt.pem; do
      if [ -f "$file" ]; then
        filename=$(basename -- "$file")
        filename_no_ext="${filename%.*}"
        new_filename="${DOCKER_CERTS_DIR}/${filename_no_ext}.crt"
        cp "$file" "$new_filename"
        echo "Copied $filename to $new_filename"
      fi
    done
    ${DOCKER_CERTS_UPDATE_COMMAND}
  fi

  echo "sccache initial state"
  ${SCCACHE} --show-stats
fi

if [ -z "${NPROC}" ]; then
  NPROC=$(nproc)
  echo "NPROC unset, setting to ${NPROC}..."
fi
echo "Using ${NPROC} cores."

if [[ "${RUN_INTEGRATED_TESTS}" = true ]]; then
  echo "Running the integrated tests has been requested."
  # We install the python environment required by ATS to run the integrated tests.
  or_die apt-get update
  or_die apt-get install -y virtualenv python3-dev python-is-python3
  ATS_PYTHON_HOME=/tmp/run_integrated_tests_virtualenv
  or_die virtualenv ${ATS_PYTHON_HOME}

  python3 -m pip cache purge

  # Setup a temporary directory to hold tests
  tempdir=$(mktemp -d)
  echo "Setting up a temporary directory to hold tests and baselines: $tempdir"
  trap "rm -rf $tempdir" EXIT
  ATS_BASELINE_DIR=$tempdir/GEOS_integratedTests_baselines
  ATS_WORKING_DIR=$tempdir/GEOS_integratedTests_working

  export ATS_FILTER="np<=32"
  ATS_CMAKE_ARGS="-DATS_ARGUMENTS=\"--machine openmpi --ats openmpi_mpirun=/usr/bin/mpirun --ats openmpi_args=--allow-run-as-root --ats openmpi_procspernode=32 --ats openmpi_maxprocs=32\" -DPython3_ROOT_DIR=${ATS_PYTHON_HOME} -DATS_BASELINE_DIR=${ATS_BASELINE_DIR} -DATS_WORKING_DIR=${ATS_WORKING_DIR}"
fi


if [[ "${CODE_COVERAGE}" = true ]]; then
  or_die apt-get update
  or_die apt-get install -y lcov
fi



# The -DBLT_MPI_COMMAND_APPEND="--allow-run-as-root;--oversubscribe" option is added for OpenMPI.
#
# OpenMPI prevents from running as `root` user by default.
# And by default user is `root` in a docker container.
# Using this option therefore offers a minimal and convenient way to run the unit tests.
#
# The option `--oversubscribe` tells OpenMPI to allow more MPI ranks than the node has cores.
# This is needed because our unit test `blt_mpi_smoke` is run in parallel with _hard coded_ 4 ranks.
# While some of our ci nodes may have less cores available.
# 
# In case we have more powerful nodes, consider removing `--oversubscribe` and use `--use-hwthread-cpus` instead.
# This will tells OpenMPI to discover the number of hardware threads on the node,
# and use that as the number of slots available. (There is a distinction between threads and cores).
GEOSX_BUILD_DIR=/tmp/geos-build
or_die python3 scripts/config-build.py \
               -hc ${HOST_CONFIG} \
               -bt ${CMAKE_BUILD_TYPE} \
               -bp ${GEOSX_BUILD_DIR} \
               -ip ${GEOSX_DIR} \
               --ninja \
               -DBLT_MPI_COMMAND_APPEND='"--allow-run-as-root;--oversubscribe"' \
               -DGEOSX_INSTALL_SCHEMA=${GEOSX_INSTALL_SCHEMA} \
               -DENABLE_COVERAGE=$([[ "${CODE_COVERAGE}" = true ]] && echo 1 || echo 0) \
               ${SCCACHE_CMAKE_ARGS} \
               ${ATS_CMAKE_ARGS}

# The configuration step is now over, we can now move to the build directory for the build!
or_die cd ${GEOSX_BUILD_DIR}

# Code style check
if [[ "${TEST_CODE_STYLE}" = true ]]; then
  or_die ctest --output-on-failure -R "testUncrustifyCheck"
  exit 0
fi

# Documentation check
if [[ "${TEST_DOCUMENTATION}" = true ]]; then
  or_die ctest --output-on-failure -R "testDoxygenCheck"
  exit 0
fi

# Performing the requested build.
if [[ "${BUILD_EXE_ONLY}" = true ]]; then
  or_die ninja -j $NPROC geosx
else
  or_die ninja -j $NPROC
  or_die ninja install

  if [[ ! -z "${DATA_BASENAME_WE}" ]]; then
    # Here we pack the installation.
    # The `--transform` parameter provides consistency between the tarball name and the unpacked folder.
    or_die tar czf ${DATA_EXCHANGE_DIR}/${DATA_BASENAME_WE}.tar.gz --directory=${GEOSX_TPL_DIR}/.. --transform "s|^./|${DATA_BASENAME_WE}/|" .
  fi
fi

if [[ ! -z "${SCCACHE_CREDS}" ]]; then
  echo "sccache post-build state"
  or_die ${SCCACHE} --show-adv-stats
fi

if [[ "${CODE_COVERAGE}" = true ]]; then
  or_die ninja coreComponents_coverage
  cp -r ${GEOSX_BUILD_DIR}/coreComponents_coverage.info.cleaned ${GEOS_SRC_DIR}/geos_coverage.info.cleaned
fi

# Run the unit tests (excluding previously ran checks).
if [[ "${RUN_UNIT_TESTS}" = true ]]; then
  if [ ${HOSTNAME} == 'streak.llnl.gov' ] || [ ${HOSTNAME} == 'streak2.llnl.gov' ]; then
    or_die ctest --output-on-failure -E "testUncrustifyCheck|testDoxygenCheck|testExternalSolvers"
  else
    or_die ctest --output-on-failure -E "testUncrustifyCheck|testDoxygenCheck"
  fi
fi

if [[ "${RUN_INTEGRATED_TESTS}" = true ]]; then
  # We split the process in two steps. First installing the environment, then running the tests.
  or_die ninja ats_environment
  
  # The tests are not run using ninja (`ninja --verbose ats_run`) because it swallows the output while all the simulations are running.
  # We directly use the script instead...
  echo "Available baselines:"
  ls -lR /tmp/geos/baselines

  echo "Running integrated tests..."
  integratedTests/geos_ats.sh --baselineCacheDirectory /tmp/geos/baselines
  tar -czf ${DATA_EXCHANGE_DIR}/test_logs_${DATA_BASENAME_WE}.tar.gz integratedTests/TestResults

  echo "Checking results..."
  bin/geos_ats_log_check integratedTests/TestResults/test_results.ini -y ${GEOS_SRC_DIR}/.integrated_tests.yaml &> $tempdir/log_check.txt
  cat $tempdir/log_check.txt

  if grep -q "Overall status: PASSED" "$tempdir/log_check.txt"; then
    echo "IntegratedTests passed. No rebaseline required."
    INTEGRATED_TEST_EXIT_STATUS=0
  else
    echo "IntegratedTests failed. Rebaseline is required."

    # Rebaseline and pack into an archive
    echo "Rebaselining..."
    integratedTests/geos_ats.sh -a rebaselinefailed

    echo "Packing baselines..."
    integratedTests/geos_ats.sh -a pack_baselines --baselineArchiveName ${DATA_EXCHANGE_DIR}/baseline_${DATA_BASENAME_WE}.tar.gz --baselineCacheDirectory /tmp/geos/baselines
    INTEGRATED_TEST_EXIT_STATUS=1
  fi

  echo "Done!"

  # INTEGRATED_TEST_EXIT_STATUS=$?
  echo "The return code of the integrated tests is ${INTEGRATED_TEST_EXIT_STATUS}"
fi

# Cleaning the build directory.
or_die ninja clean


# Clean the repository
or_die cd ${GEOS_SRC_DIR}/inputFiles
find . -name *.pyc | xargs rm -f


# If we're here, either everything went OK or we have to deal with the integrated tests manually.
if [[ ! -z "${INTEGRATED_TEST_EXIT_STATUS+x}" ]]; then
  echo "Exiting the build process with exit status ${INTEGRATED_TEST_EXIT_STATUS} from the integrated tests."
  exit ${INTEGRATED_TEST_EXIT_STATUS}
else
  echo "Exiting the build process with exit status 0."
  exit 0
fi