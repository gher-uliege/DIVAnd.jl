# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot

include("./prep_dirs.jl")

# observations
x = [2  3  4];
y = [2  3  4];
f=ones(3)

# final grid
xi,yi = ndgrid(range(0,stop=6,length=50),range(0,stop=5,length=30));


# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);

# correlation length
len = 1;

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field
fireg,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);

@show sampler1=DIVAnd_sampler((pm,pn),len)

fis,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,(len*0.5,len*1.5),epsilon2);

@show  sampler1=DIVAnd_sampler((pm,pn),(len*0.5,len*1.5))

pm = ones(size(xi)) ./ ((1 .+ xi/5).*(xi[2,1]-xi[1,1]));
pn = ones(size(yi)) ./ ((1 .+ yi/5).*(yi[1,2]-yi[1,1]));

lx = 0.5 .+ xi/5
ly = 0.5 .+ yi/5

finu,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,(ly,lx),epsilon2);
@show  sampler1=DIVAnd_sampler((pm,pn),(lx,ly))

# plotting of results
subplot(1,3,1);
pcolor(xi,yi,fireg);
colorbar()
clim(0,1)
plot(x,y,"k.");
subplot(1,3,2);
pcolor(xi,yi,fis);
colorbar()
clim(0,1)
subplot(1,3,3);
pcolor(xi,yi,finu);
colorbar()
clim(0,1)

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => ".png")));
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
