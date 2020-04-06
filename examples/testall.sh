for i in *.jl; do sbatch  $(awk  /#SBATCH/ { print } $i) --job-name=${i%.jl} --output=${i%.jl}-%j.log submit_julia_1cpu.sh $i; done
