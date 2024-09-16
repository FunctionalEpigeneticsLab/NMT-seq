#!/usr/bin/bash
#SBATCH -t 03:00:00
#SBATCH --mem=10G
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

/usr/bin/singularity build nmt_rna_pipeline.sif nmt_rna_pipeline.def
