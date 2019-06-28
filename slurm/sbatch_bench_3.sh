#!/bin/bash

#SBATCH -J "3_benchmark"
#SBATCH --array=0-46999%4096
#SBATCH --account=p_mwrc
#SBATCH --time=1:59:59
#SBATCH --mem-per-cpu=2583
#SBATCH --partition=haswell
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/p_mwrc/log/%x-wp_%a-job_%A.out

export JOB_BASEDIR="/scratch/p_mwrc/${SLURM_JOB_NAME}"
export JOB_HPC_SAVEDIR="${JOB_BASEDIR}/results"

if [ ! -d ${JOB_BASEDIR} ]; then
	mkdir ${JOB_BASEDIR}
	mkdir ${JOB_HPC_SAVEDIR}
fi

NC=1000
JOBS=('4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '34' '35' '36' '37' '38' '39' '40' '41' '42' '43' '44' '45' '46' '47' '48' '49' '50')
CIDX=$((${SLURM_ARRAY_TASK_ID} % ${NC}))
JID=$((${SLURM_ARRAY_TASK_ID} / ${NC}))

EXE=$(readlink -e ../code/bench_${JOBS[$JID]}_3)
DATA=$(readlink -e ../data/wp-bench.h5)

module load HDF5/1.10.5-foss-2019a
module load foss
module load GCC
module load mosek
module load Gurobi

# for better error reporting (ticket 2018071041004097)
echo $LD_LIBRARY_PATH
ls -la $EBROOTGUROBI/lib/libgurobi80.so

srun $EXE $DATA $CIDX
