#!/usr/bin/bash
#SBATCH -t 00:20:00
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cluster wice
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -A lp_big_wice_cpu
#SBATCH --partition dedicated_big_bigmem


export SINGULARITY_TMPDIR=$VSC_SCRATCH/singularity_tmp
mkdir -p $SINGULARITY_TMPDIR
export SINGULARITY_CACHEDIR=$VSC_SCRATCH/singularity_cache
mkdir -p $SINGULARITY_CACHEDIR

singularity build allc_env.sif allc_env.def
