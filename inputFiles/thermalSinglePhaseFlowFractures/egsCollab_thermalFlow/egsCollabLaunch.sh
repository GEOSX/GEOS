#!/bin/tcsh
##### These lines are for Slurm
#SBATCH -N 1
#SBATCH -J egsCollab
#SBATCH -t 2:00:00
#SBATCH -p pbatch
#SBATCH --mail-type=ALL
#SBATCH -A vortex
#SBATCH -o /usr/workspace/cusini1/geosx/geosx_integratedTests/GEOSX/egsCollab.out

##### These are shell commands
date
cd /usr/workspace/cusini1/geosx/geosx_integratedTests/GEOSX/build-quartz-gcc@12-release
##### Launch parallel job using srun
setenv OMP_NUM_THREADS 1
srun -n8 /bin/geosx -i ../inputFiles/thermalSinglePhaseFlowFractures/egsCollab_thermalFlow/egsCollab_thermalFlow_initialCond_fine.xml -x 2 -y 2 -z 2 -o Output_collab
echo 'Done'