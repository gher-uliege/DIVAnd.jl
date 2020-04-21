# test all tests with SLURM starting with DIVAnd

for i in DIVAnd*.jl; do
    echo submit $i
    sbatch  $(awk  '/#SBATCH/ { print $2 }'  $i) --job-name=${i%.jl} --output=${i%.jl}-%j.log submit_julia_1cpu.sh $i; 
done
