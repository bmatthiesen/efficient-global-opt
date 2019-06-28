#!/bin/sh

#SBATCH -J "afsnd"
#SBATCH --array=790,1330,2125,2500,2680,2725,2905,3655,3745,4540,4630,4720,5245,5395,5515,5560,5665,6025,6100,6325,6910,7000,7030,7315,7525,7900,8725,8770,8875,8995,9040,9160,9505,9640,9820,9985,10105,10150,10390,10480,10585,10660,10750,10810,10840,10870,11020,11035,11440,11530,11815,11830,11995,12055,12700,12715,12850,12910,13360,13375,13465,13675,13780,13975,14140,14170,14215,14290,14365,14695,14995
#SBATCH --account=p_mwrc
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=10583
#SBATCH --partition=haswell
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/p_mwrc/log/%x-wp_%a-job_%A.out

export JOB_BASEDIR="/scratch/p_mwrc/${SLURM_JOB_NAME}"
export JOB_HPC_SAVEDIR="${JOB_BASEDIR}/results"
export JOB_CKPTDIR="${JOB_BASEDIR}/ckpt"

if [ ! -d ${JOB_BASEDIR} ]; then
	mkdir ${JOB_BASEDIR}
	mkdir ${JOB_HPC_SAVEDIR}
	mkdir ${JOB_CKPTDIR}
fi

export JOB_WD="${JOB_CKPTDIR}/wp${SLURM_ARRAY_TASK_ID}"

if [ ! -d ${JOB_WD} ]; then
	mkdir ${JOB_WD}
fi

EXE=$(readlink -e ../code/afsnd2)
DATA=$(readlink -e ../data/wp.h5)
INTERVAL=14400

cd ${JOB_WD}

module load foss
module load GCC
module load mosek
module load Gurobi
module load DMTCP

source ${DMTCP_ROOT}/bin/bash
start_coordinator -i $INTERVAL

# for better error reporting (ticket 2018071041004097)
echo $LD_LIBRARY_PATH
ls -la $EBROOTGUROBI/lib/libgurobi80.so

dmtcp_launch $EXE $DATA
