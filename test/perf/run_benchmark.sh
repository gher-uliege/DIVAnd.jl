#!/bin/bash
#SBATCH --job-name=julia
#SBATCH --output=output-%N-%j.out
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50000

# run as:
# sbatch run_benchmark.sh benchmark1.jl 4 varanalysis
# for i in $(seq 2 4); do sbatch run_benchmark.sh benchmark1.jl $i varanalysis; done

module load  EasyBuild  Python/3.5.1-foss-2016a
srun julia "$@"
