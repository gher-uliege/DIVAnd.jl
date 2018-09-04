#SBATCH --mem-per-cpu=8000


# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot

if VERSION >= v"0.7"
    using LinearAlgebra
    using Random
    using Statistics
end

include("./prep_dirs.jl")

if VERSION >= v"0.7"
   Random.seed!(1234)
else
   srand(1234)
end
# observations
nobs=100
x = 0.01 .+ 0.98*rand(nobs);
y = 0.01 .+ 0.98*rand(nobs);
x = -0.1 .+ 1.2*rand(nobs);
y = -0.1 .+ 1.2*rand(nobs);
f = sin.(x*6) .* cos.(y*6);
#f=-1+2*x
# x=[0.5,0.75]
# y=[0.5,0.75]
# f=[1,1]

# final grid
xi,yi = ndgrid(range(0,stop=1,length=950),range(0,stop=1,length=830));

# reference field
fref = sin.(xi*6) .* cos.(yi*6);

# all points are valid points
mask = trues(size(xi));

mask[300:600,400:600] .= false

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);

# correlation length
len = 0.03;

# obs. error variance normalized by the background error variance
epsilon2 = 1.;
vscale=0.001
vscale=0

# fi is the interpolated field
@time fiexOLD,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,(len,0.5*len),epsilon2;alphabc=0);

residue=DIVAnd_residualobs(s,fiexOLD)

@time fiOLD,errOLD,residueGO = DIVAndgo(mask,(pm,pn),(xi,yi),(x,y),f,(len,0.5*len),epsilon2;alphabc=0);


@show size(residue),size(residueGO)
@show norm(residue-residueGO)

figure("P1")

subplot(1,3,1)
title("Old version")
pcolor(xi,yi,fiOLD)
clim(-1,1)
subplot(1,3,2)
title("Reference")
pcolor(xi,yi,fiexOLD)
clim(-1,1)
subplot(1,3,3)
rms=sqrt(var(fiexOLD-fiOLD))
title("rms $rms")
pcolor(xi,yi,fiexOLD-fiOLD)
clim(-0.05,0.05)
plot(x,y,"k.");
colorbar()

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_1.png")));
savefig(figname)
@info "Saved figure as " * figname

@time fiex,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,(len,0.5*len),epsilon2;alphabc=1.);

residue=DIVAnd_residualobs(s,fiex)

@time fi,erri,residueGO = DIVAndgo(mask,(pm,pn),(xi,yi),(x,y),f,(len,0.5*len),epsilon2;alphabc=1.);

@show size(residue),size(residueGO)

@show norm(residue-residueGO)

figure("Pp")

subplot(1,3,1)
title("NEW version")
pcolor(xi,yi,fi)
clim(-1,1)
subplot(1,3,2)
title("Reference")
pcolor(xi,yi,fiex)
clim(-1,1)
subplot(1,3,3)
rms=sqrt(var(fiex-fi))
title("rms $rms")
pcolor(xi,yi,fiex-fi)
clim(-0.05,0.05)
plot(x,y,"k.");
colorbar()

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_2.png")));
savefig(figname)
@info "Saved figure as " * figname

figure("Pppp")

subplot(1,2,1)
title("cpme NEW version")
pcolor(xi,yi,erri)
clim(0,1)
subplot(1,2,2)
title("cpme OLD version")
pcolor(xi,yi,errOLD)
clim(0,1)
colorbar()

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_3.png")));
savefig(figname)
@info "Saved figure as " * figname

figure()
scatter(residue,residueGO)

# Copyright (C) 2014, 2018 Alexander Barth         <a.barth@ulg.ac.be>
#                          Jean-Marie Beckers   <JM.Beckers@ulg.ac.be>
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
