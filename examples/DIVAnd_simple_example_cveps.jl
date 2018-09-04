# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot
if VERSION >= v"0.7"
    using Random
end

include("./prep_dirs.jl")

if VERSION >= v"0.7"
   Random.seed!(1234)
else
   srand(1234)
end
# observations
nobs = 90
x = rand(nobs);
y = rand(nobs);
f = sin.(x*6) .* cos.(y*6);
f = f + randn(nobs);

# final grid
xi,yi = ndgrid(range(0,stop=1,length=70),range(0,stop=1,length=70));

# reference field
fref = sin.(6xi) .* cos.(6yi);

# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);

# correlation length
len = 0.2;

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field
bestfactore = 1
cvval = 99999

figure()

for imeth=0:3
    bestfactore, cvval,cvvalues, x2Ddata,cvinter,xi2D = DIVAnd_cv(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2,0,4,imeth);

    subplot(2,2,imeth+1)
    plot(xi2D,cvinter,"-")
    xlabel("Log10 scale factor e2")
    plot(x2Ddata,cvvalues,".")
    plot(log10.(bestfactore), cvval,"o")
    title("Method $imeth")
end

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_1.png")))
savefig(figname)
@info "Created figure " * figname

# De Rosier type of approach

figure("another one")

epsbest1=epsilon2*bestfactore
cvbest1=cvval

cvbest2=zeros(20);
eps2=zeros(20)

for i=1:20
    global epsilon2
    cvval,factor=DIVAnd_cv(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2,0,0,3);
    eps2[i]=epsilon2;
    cvbest2[i]=cvval;
    epsilon2=epsilon2*factor
end

epsilon2
epsbest1
cvbest1
cvbest2
eps2

plot(log10.(eps2),cvbest2,".",log10.(epsbest1),cvbest1,"o",log10.(eps2[end]),cvbest2[end],"+")

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_2.png")))
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
