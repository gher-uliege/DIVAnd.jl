# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range, argmin
using PyPlot
if VERSION >= v"0.7"
    using LinearAlgebra
    using Statistics
end


include("./prep_dirs.jl")

# observations
x = [-0.7,0.7,0.7,-0.7];
y = [0.,-.7,0.7,0.7];
f = [1.,-1.,1.,-1.];

# final grid
xi,yi = ndgrid(range(-1,stop=1,length=81),range(-1,stop=1,length=81));
xifin,yifin = ndgrid(range(-10,stop=10,length=801),range(-10,stop=10,length=801));

varb1=0
Bi=0
Bold=0

pmfin = ones(size(xifin)) / (xifin[2,1]-xifin[1,1]);
pnfin = ones(size(xifin)) / (yifin[1,2]-yifin[1,1]);

# correlation length
len = 1/4.1*8/9.75;
@show len
@show 2/len

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

maskfin = trues(size(xifin));
# fi is the interpolated field
fifin,sfin = DIVAndrun(maskfin,(pmfin,pnfin),(xifin,yifin),(x,y),f,len,epsilon2;alphabc=0);

firef=fifin[401-40:401+40,401-40:401+40];


# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);

@show pm[1,1]*len
# correlation length


# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field
fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;alphabc=1.09);

Bnew=diag(s.P)
@show mean(Bnew)

fiold,sold = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,1E6;alphabc=0);
Bold=diag(sold.P);
@show mean(Bold)

#pcolor(reshape(Bold,81,81));colorbar()
figure("ww")

fiold,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;alphabc=0);

@show sqrt(var(fiold-fi))/sqrt(var(fiold-firef))
varr=zeros(100)
rms=zeros(100)
al=zeros(100)
for ii=1:100
    alen=0.25+ii/100*1.25
    fi,si = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;alphabc=alen);
    fbi,sbi = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,1E6;alphabc=alen);
    Bi=diag(sbi.P);
    varr[ii]=sqrt(var(Bi))/sqrt(var(Bold))
    rms[ii]=sqrt(var(firef-fi))/sqrt(var(fiold-firef))
    al[ii]=alen

end


figure("min")
subplot(1,2,1)
plot(al,varr,"-")
subplot(1,2,2)
plot(al,rms,"-")
@show al[argmin(varr)]

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_1.png")));
savefig(figname)
@info "Saved figure as " * figname


alphaopt=al[argmin(varr)]
fi,swhat = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;alphabc=alphaopt);
#fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;alphabc=0);
fidir,swhat = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;alphabc=2.06);

figure("other")

# plotting of results
subplot(2,2,1);
rmval=sqrt(var(firef-firef))
title("Reference rms $rmval")
pcolor(xi,yi,firef);
colorbar()
clim(-1,1)
plot(x,y,"k.");

subplot(2,2,2);
rmval=sqrt(var(fi-firef))
title("optimal alpha  rms $rmval")
pcolor(xi,yi,fi);
colorbar()
clim(-1,1)

subplot(2,2,3);
rmval=sqrt(var(fidir-firef))
title("alpha=1.09  rms $rmval")
pcolor(xi,yi,fidir);
colorbar()
clim(-1,1)

subplot(2,2,4);
rmval=sqrt(var(fiold-firef))
title("old bc  rms $rmval")
pcolor(xi,yi,fiold);
colorbar()
clim(-1,1)

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_2.png")));
savefig(figname)
@info "Saved figure as " * figname


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
