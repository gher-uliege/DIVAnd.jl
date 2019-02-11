using DIVAnd
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

sz = (80,90,3)
f = fill(NaN,sz);
f[24:40,30:65,1] .= 1;
f[30:50,60:75,2] .= 2;
f[20:35,20:55,3] .= 1.5;

L = 4
f2 = DIVAnd.hmerge(f,L);

@test all(1 .<= f2[isfinite.(f2)] .<= 2)
#pcolor(copy(f2'));
