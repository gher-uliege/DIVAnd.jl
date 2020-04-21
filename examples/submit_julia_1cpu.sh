#!/bin/bash
#SBATCH --output=output-%N-%j.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000

export script="$1"
echo Running script $script

bt0=$(date +%s)

unset DISPLAY

julia "$script"

bt1=$(date +%s)

awk  " BEGIN { print \"Run time (hours): \",($bt1 - $bt0)/3600 } "



