#!/bin/bash
#SBATCH --job-name=julia
#SBATCH --mail-type=ALL
#SBATCH --output=output-%N-%j.out
#SBATCH --time=10:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=18000

# run as run_benchmark.sh benchmark.jl

module load  EasyBuild  Python/3.5.1-foss-2016a
srun --mem-per-cpu=16000 --ntasks=1 --time=10:20:00 julia "$@"
