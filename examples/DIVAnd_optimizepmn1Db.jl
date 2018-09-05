#
using DIVAnd
using Compat: @info, range, argmin
using PyPlot
if VERSION >= v"0.7"
    using LinearAlgebra
    using Statistics
end

include("./prep_dirs.jl")

xi=0
xiref=0

aj=zeros(300)
vj=zeros(300)

# Calculate reference solution on a very wide domain

xiref = collect(range(-100.0,stop=100.0,length=2001));
len=2.
x = [-3.0, 8.];
f = [1.,  1.];
pmref = ones(size(xiref)) / (xiref[2]-xiref[1]);
maskref = trues(size(xiref));
epsilon2=1.
firef,sref = DIVAndrun(maskref,(pmref,),(xiref,),(x,),f,len,epsilon2,alphabc=0);

epsilon2large = 10000.;
firefb,s = DIVAndrun(maskref,(pmref,),(xiref,),(x,),f,len,epsilon2large,alphabc=0);
bref=diag(s.P)

xirefl = collect(range(-10.0,stop=10.0,length=201));
len=2
x = [-3.0, 8.];
f = [1.,  1.];
pmrefl = ones(size(xirefl)) / (xirefl[2]-xirefl[1]);
maskrefl = trues(size(xirefl));
epsilon2=1.
firefl,srefl = DIVAndrun(maskrefl,(pmrefl,),(xirefl,),(x,),f,len,epsilon2,alphabc=0);

firefbb,s = DIVAndrun(maskrefl,(pmrefl,),(xirefl,),(x,),f,len,epsilon2large,alphabc=0);
brefl=diag(s.P)

figure("Reference")
rmsdiff=sqrt(var(firef[901:1101]-firefl))
title("Solution in infinite domain and finite domain, rms = $rmsdiff")
plot(xiref[801:1201],firef[801:1201],"-",xirefl,firefl,".")
figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_reference.png")))
savefig(figname);
@info "Created figure " * figname
# Now try to optimize BC
aj=zeros(500)
vj=zeros(500)
rj=zeros(500)

xi = collect(range(-10.0,stop=10.0,length=201));

# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2]-xi[1]);
# obs. error variance normalized by the background error variance



for j=1:500

    alen=0.25+j/100

    #Test to push boundary to wider distance:

    # pm[201]=1./(alen*len);
    # pm[1]  =1./(alen*len);


    # fi is the interpolated field
    fi2,s = DIVAndrun(mask,(pm,),(xi,),(x,),f,len,epsilon2large,alphabc=alen);




    aj[j]=alen
    vj[j]=var(diag(s.P))

    # Now with real data for comparison of analysis
    fi2,s = DIVAndrun(mask,(pm,),(xi,),(x,),f,len,epsilon2,alphabc=alen);

    rj[j]=sqrt(var(firef[901:1101]-fi2))

end


@show argmin(vj)
@show argmin(rj)

alpha=aj[argmin(vj)]
varb=vj[argmin(vj)]

figure("Optimization")

subplot(1,2,1)
title("variance of diag(B)")
plot(aj,vj,"-")
subplot(1,2,2)
title("rms(reference-analysis)")
plot(aj,rj,"-")
figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_optimisation.png")))
savefig(figname);
@info "Created figure " * figname

# Finally solution with optimized parameter

figure("Optimal solution 1")
alen=aj[argmin(vj)]
@show alen


# pm[201]=1./(alen*len);
# pm[1]  =1./(alen*len);
# Now with real data for comparison of analysis
fi2,s = DIVAndrun(mask,(pm,),(xi,),(x,),f,len,epsilon2,alphabc=alen);
rmsdiff=sqrt(var(firef[901:1101]-fi2))
fi2b,s = DIVAndrun(mask,(pm,),(xi,),(x,),f,len,epsilon2large,alphabc=alen);

subplot(2,1,1)
title("Solution in infinite domain and modified finite domain, Bversion , rms = $rmsdiff")
plot(xiref[801:1201],firef[801:1201],"-",xi,fi2,".")

subplot(2,1,2)
bi=diag(s.P)
title("B in infinite domain, finite domain and modified finite domain, Bversion")
plot(xiref[801:1201],bref[801:1201],"-",xi,bi,".",xi,brefl,".")
figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_optimal1.png")))
savefig(figname);
@info "Created figure " * figname

figure("Optimal solution 2")
alen=aj[argmin(rj)]
@show alen


# pm[201]=1./(alen*len);
# pm[1]  =1./(alen*len);
# Now with real data for comparison of analysis
fi2,s = DIVAndrun(mask,(pm,),(xi,),(x,),f,len,epsilon2,alphabc=alen);
rmsdiff=sqrt(var(firef[901:1101]-fi2))
fi2b,s = DIVAndrun(mask,(pm,),(xi,),(x,),f,len,epsilon2large,alphabc=alen);
subplot(2,1,1)
title("Solution in infinite domain and modified finite domain, rmsversion, rms = $rmsdiff", fontsize=14)
plot(xiref[801:1201],firef[801:1201],"-",xi,fi2,".")
subplot(2,1,2)

bi=diag(s.P)
title("B in infinite domain, finite domain and modified finite domain, rmsversion", fontsize=14)
plot(xiref[801:1201],bref[801:1201],"-",xi,bi,".",xi,brefl,".")
figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_optimal2.png")))
savefig(figname)
@info "Created figure " * figname

# Copyright (C) 2014, 2018 Alexander Barth <a.barth@ulg.ac.be>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, see <http://www.gnu.org/licenses/>.
