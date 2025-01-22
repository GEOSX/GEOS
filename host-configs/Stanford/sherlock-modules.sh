module --force purge 
module load devel
module load math
module load system
module load git
module load git-lfs
module load cmake
module load gcc/10.1.0
module load openmpi/4.1.2
module load openblas/0.3.10
module load libevent/2.1.12 # needed for silo
module load py-mpi4py/3.1.3_py39
module load py-h5py/3.7.0_py39
module unload cuda
module load cuda/12.4.0 
module load ninja/1.9.0

