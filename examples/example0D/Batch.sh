#!/bin/bash
#SBATCH --job-name="qbox-ucriverside"
#SBATCH --output="riverside.out"
#SBATCH --nodes=1
#SBATCH --account=m2956
#SBATCH -t 02:59:00
#SBATCH -q regular
#SBATCH -C gpu



#module swap PrgEnv-gnu PrgEnv-nvidia #PrgEnv-nvhpc
module load PrgEnv-nvidia
module load cray-fftw cray-libsci
module unload darshan

export MPICH_GPU_SUPPORT_ENABLED=0

PROFILE=0

#export QB_HOME=/global/cfs/cdirs/m2956/adrianpd/Lastest_dev/src
export QB_HOME=/global/homes/s/sliu424/QBOX_TDDFT_GPU/src


if [[ $PROFILE -eq 1 ]]
then
        module load perftools
        QB=$QB_HOME/qb+pat
else
        QB=$QB_HOME/qb
fi



#export XERCES_HOME=/global/cfs/cdirs/m2956/adrianpd/xerces-perl
export XERCES_HOME=/global/homes/s/sliu424/xerces-3.2.4
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$XERCES_HOME/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$XERCES_HOME/lib


export OMP_NUM_THREADS=1

ppn=$((SLURM_CPUS_ON_NODE/OMP_NUM_THREADS))
srun -n 8 -c $OMP_NUM_THREADS --gpus-per-node=4 $QB -nstb 8 -nkpb 1 -nspb 1 < td.in > output_vector.out




#srun -n 8 -c $OMP_NUM_THREADS --gpus-per-node=4 nsys profile -t nvtx,cuda --stats=true --force-overwrite true  $QB -nstb 8 -nkpb 1 -nspb 1 < td.in > salida_gpu_profile.out


