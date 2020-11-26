using DIVAnd
using Random
using Test

# Observations
np = 20
nt = 5

obslon = rand(np);
obslat = rand(np);
obsval = rand(np, nt);

# Grid
lonr = 0.:0.05:1.
latr = 0.:0.1:1.
mask, (pm, pn), (xi, yi) = DIVAnd_rectdom(lonr, latr);

# Parameters
len = 0.2
epsilon2 = 5.

# First analysis
@time fi1ref,s = DIVAndrun((mask),(pm,pn),
    (xi,yi),(obslon,obslat),Float64.(obsval[:,1]),len,epsilon2);

# Other analysis using structure 's'
fi = Array{Float64, 3}(undef, length(lonr), length(latr), nt)
for i = 1:nt
    fpi = s.P * (s.H' * (s.R \ obsval[:,i]))
    fi[:,:,i] = unpack(s.sv, fpi, NaN)[1];
end

# Test with the first time step
@test fi1ref == fi[:,:,1]
