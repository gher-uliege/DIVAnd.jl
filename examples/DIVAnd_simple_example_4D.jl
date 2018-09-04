# A simple example of DIVAnd in 4 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot

include("./prep_dirs.jl")

# function to interpolate
fun(x,y,z,t) = sin.(6x) .* cos.(6y)+sin.(6z) .* cos.(6x) .* sin.(2*pi*t) ;

# observations
nobs = 200;
x = rand(nobs);
y = rand(nobs);
z = 0.5 .+ 0.01*rand(nobs);
t = rand(nobs);
f = fun.(x,y,z,t)

# final grid
testsizexy=100
testsizez=2
testsizet=4

# mask: all points are valid points
# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension,...
mask,(pm,pn,po,pq),(xi,yi,zi,ti) = DIVAnd_rectdom(range(0,stop=1,length=testsizexy),
                                                  range(0,stop=1,length=testsizexy),
                                                  range(0,stop=1,length=testsizez),
                                                  range(0,stop=1,length=testsizet))

# reference field
fref = fun.(xi,yi,zi,ti)

# correlation length
len = (0.1,0.1,0.1,0.1);

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field
@time fi,s = DIVAndrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2; moddim=[0,0,0,1]);

# plotting of results
subplot(1,2,1);
pcolor(xi[:,:,1,3],yi[:,:,1,3],fref[:,:,1,3],vmin=-1.,vmax=1.);
clim(-1,1)
plot(x,y,"k.");
title("Analytical field and obs.");

subplot(1,2,2);
pcolor(xi[:,:,1,3],yi[:,:,1,3],fi[:,:,1,3],vmin=-1.,vmax=1.);
colorbar()
clim(-1,1)
title("Interpolated field");

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => ".png")))
savefig(figname)
@info "Created figure " * figname

# Copyright (C) 2014, 2018 Alexander Barth <a.barth@ulg.ac.be>
#                          Jean-Marie Beckers <JM.Beckers@ulg.ac.be>
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
