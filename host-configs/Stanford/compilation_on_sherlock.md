# GEOS Compilation on Sherlock

## Overview
This guide provides a step-by-step process for compiling GEOS simulator on the Stanford Sherlock cluster. The compilation involves both Third-Party Libraries (TPLs) and GEOS. 

### Important Note

The shell scripts `clone.sh`, `tpls.sh`, `geos.sh` (partial scripts) and `compile_geos.sh` (main script) must be created and populated with the relevant content outlined in this document.

To facilitate this, please follow these organizational instructions:

* In your working directory, create a folder named `build_utils`;
* Place the partial scripts  `clone.sh`, `tpls.sh`, `geos.sh` within the build_utils directory;
* Ensure that the main script `compile_geos.sh` is located in your working directory at the same level as the `build_utils` folder.

This organizational structure is relevant to ensure the successful execution of the compilation process. 

The following illustrates the files structure. If you run `ls` and `ls -l build_utils/` commands in your working directory, you should see something like this:

```
[suid@sh04-ln04 login /home/groups/pi_suid/suid]$ ls
build_utils  compile_geos.sh

[suid@sh04-ln04 login /home/groups/pi_suid/suid]$ ls -l build_utils/
total 96
-rw-r--r-- 1 suid pi_suid  565 Jan 24 18:07 clone.sh
-rw-r--r-- 1 suid pi_suid 1424 Jan 24 18:46 geos.sh
-rw-r--r-- 1 suid pi_suid 1380 Jan 24 18:46 tpls.sh
```

## Compilation Steps

### Main Step: Execute the Compile Script
The compile_geos.sh script executes the compilation process:

```bash
source compile_geos.sh
```

This script orchestrates the compilation process by calling other scripts in a defined order:
1. Clone the Sources
2. Compile TPLs
3. Compile GEOS

### Step 1: Clone the Sources
The **`clone.sh`** script is responsible for fetching the required repositories and initializing the necessary submodules. The script performs the following actions:

- It loads the appropriate modules for Git.
- It clones the Third-Party Libraries (`thirdPartyLibs`) and the GEOS source code.
- It initializes the Git Large File Storage (LFS) and updates the submodules.

**Content of `clone.sh`:**

```bash
#!/bin/bash
# Load necessary modules
module load system
module load git/2.45.1 git-lfs/2.4.0

# Clone the sources
GIT_CLONE_PROTECTION_ACTIVE=false git clone https://github.com/GEOS-DEV/thirdPartyLibs.git
cd thirdPartyLibs || { echo "Failed to enter thirdPartyLibs directory"; exit 1; }
git lfs install
git pull
git submodule init
git submodule update
cd ..

GIT_CLONE_PROTECTION_ACTIVE=false git clone https://github.com/GEOS-DEV/GEOS.git
cd GEOS || { echo "Failed to enter GEOS directory"; exit 1; }
git lfs install
git submodule init
git submodule update
cd ..
```

### Step 2: Compile TPLs
The **`tpls.sh`** script configures and compiles the Third-Party Libraries. This involves copying a custom configuration file and executing the build commands:

- It loads the necessary modules.
- It copies the CMake configuration file (`sherlock-custom.cmake`) into the appropriate directory for GEOS.
- It executes the `config-build.py` script to configure TPLs for Debug builds before running `make` to compile them.

**Content of `tpls.sh`:**

```bash
#!/bin/bash
#SBATCH --job-name=tpls_job        # Name of the job
#SBATCH --output=tpls_output.log  # Output log file 
#SBATCH --error=tpls_error.log    # Error log file 
#SBATCH --nodes=1                  # Use one node
#SBATCH --ntasks=1                 # Number of tasks (usually for MPI, set to 1 for non-MPI)
#SBATCH --cpus-per-task=4          # Request 4 CPU cores
#SBATCH --mem=16G                   # Request 16 GB of memory
#SBATCH --time=02:00:00            # Set a time limit of 2.0 hours
#SBATCH --partition=dev            # Specify the partition

# Email notifications
#SBATCH --mail-type=END,FAIL       # Email notifications for job completion and failure
#SBATCH --mail-user=suid@stanford.edu    # Replace with your email address

# Load necessary modules
source GEOS/host-configs/Stanford/sherlock-modules.sh

# Step 2: Configure TPLs
cd thirdPartyLibs/ || { echo "Failed to enter thirdPartyLibs directory"; exit 1; }
python3 scripts/config-build.py -hc ../GEOS/host-configs/Stanford/sherlock-gcc10.cmake -bt Debug -DNUM_PROC=4

# Step 3: Compile TPLs Debug
cd build-sherlock-gcc10-debug/ || { echo "Failed to enter build-sherlock-gcc10-debug directory"; exit 1; }
make
cd ../..
```

The above is an employing a CMake configuration file `sherlock-gcc10.cmake` wich is test in the continuous integratoin pipline. This file maps some of the loaded modules to configure TPLs (Third-Party Libraries) and GEOS.

The aforementioned process utilizes a CMake configuration file, sherlock-gcc10.cmake, which is tested within the continuous integration pipeline. The file `sherlock-modules.sh` load modules employed in the cmake configuration file for both TPL's and GEOS.


### Step 3: Configure GEOS
The **`geos.sh`** script takes care of configuring and compiling the GEOS simulator itself. It performs these tasks:

- It loads the necessary modules (the same as for the TPLs).
- It retrieves the absolute path for the TPLs installation.
- It executes the `config-build.py` script for GEOS, linking it to the previously built TPLs.
- Finally, it compiles GEOS using `make -j`.

**Content of `geos.sh`:**

```bash
#!/bin/bash
#SBATCH --job-name=geos_job         # Name of the job
#SBATCH --output=geos_output.log  # Output log file 
#SBATCH --error=geos_error.log    # Error log file 
#SBATCH --nodes=1                  # Use one node
#SBATCH --ntasks=1                 # Number of tasks (usually for MPI, set to 1 for non-MPI)
#SBATCH --cpus-per-task=4          # Request 4 CPU cores
#SBATCH --mem=16G                   # Request 16 GB of memory
#SBATCH --time=02:00:00            # Set a time limit of 2.0 hours
#SBATCH --partition=dev            # Specify the partition

# Email notifications
#SBATCH --mail-type=END,FAIL       # Email notifications for job completion and failure
#SBATCH --mail-user=suid@stanford.edu    # Replace with your email address

# Load necessary modules
source GEOS/host-configs/Stanford/sherlock-modules.sh

# Step 4: Configure GEOS
cd GEOS/ || { echo "Failed to enter GEOS directory"; exit 1; }

# Get absolute path for TPLs installation
tpls_path=$(realpath ../thirdPartyLibs/install-sherlock-gcc10-debug/)
python3 scripts/config-build.py -hc host-configs/Stanford/sherlock-gcc10.cmake -bt Debug -D GEOS_TPL_DIR="$tpls_path"

# Step 5: Compile GEOS Debug
cd build-sherlock-gcc10-debug/ || { echo "Failed to enter build-sherlock-gcc10-debug directory"; exit 1; }
make -j 4
cd ../..
```

## Compiling GEOS

The `compile_geos.sh` file combines the above steps into a unified process. It handles the sequence and dependencies between the jobs via the flag `--dependency`:

**Content of `compile_geos.sh`:**

```bash
# Clone sources
source build_utils/clone.sh

# Submit the first job for TPLs compilation
tpls_id=$(sbatch build_utils/tpls.sh | awk '{print $4}')

# Submit the GEOS compilation job with a dependency on the TPLs job
sbatch --dependency=afterok:$tpls_id build_utils/geos.sh
```

The GEOS compilation will be submitted only if the TPL job succeeds. This allows us to rapidly allocate small resources in the form of a `dev` partition. See [Sherlock documentation](https://www.sherlock.stanford.edu/docs/user-guide/running-jobs/?h=sh_part#available-resources) for available resources and types of partitions.

### Execution
To start the compilation process, simply run:

```bash
source compile_geos.sh
```

This will initiate the compilation while managing job dependencies. You will receive email notifications regarding job completion or failure.

### Monitoring Progress
To monitor the output while the compilation is in progress, use the following command:

For step 2:

```bash
tail -f tpls_output.log
```

For step 3:

```bash
tail -f geos_output.log
```

## Summary
This guide document the compilation process of GEOS on Stanford's Sherlock cluster. The process effectively employs SLURM's resource management functionalities to streamline jobs execution in sequence. For advanced usage or specific needs, additional configurations and modifications may be required.

## References
- [GEOSX Documentation](https://geosx-geosx.readthedocs-hosted.com/en/latest/#)
- [Sherlock Documentation](https://www.sherlock.stanford.edu/docs/)