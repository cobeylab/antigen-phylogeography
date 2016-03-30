#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=stdout.txt
#SBATCH --error=stderr.txt
#SBATCH --mem=10240
#SBATCH --time=02:00:00

let START=$SLURM_ARRAY_TASK_ID*$N_PER_JOB
let END=$START+${N_PER_JOB}-1

# if end>total tasks then end

for SUB_JOB_NUM in `seq $START $END`
do
	cd $SUB_JOB_NUM
	${ANTIGEN_ROOT}/run.py > "out.bashscreen"
	excessDiversityFile="out.tmrcaLimit"
	extinctFile="out.extinct"
	excessDiversity=${ANTIGEN_ROOT}/example_sweep/results/$SUB_JOB_NUM/$excessDiversityFile
	extinct=${ANTIGEN_ROOT}/example_sweep/results/$SUB_JOB_NUM/$extinctFile
	cd ${ANTIGEN_ROOT}/example_sweep
	if [ ! -f "$excessDiversity" -a ! -f "$extinct" ]
	then
	MathKernel -noprompt -run < antigen-analysis-noFig.m ${ANTIGEN_ROOT}/example_sweep/results/$SUB_JOB_NUM --display `hostname`:1 > "out.math"
	fi
	cd ${ANTIGEN_ROOT}/example_sweep/results/
done

# Use this instead to restart runs that go extinct:
# ${ANTIGEN_ROOT}/run_py --with-restarts
