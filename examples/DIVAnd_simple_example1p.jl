# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.


# To test a single point
#Analysis in very large domain,
#diagnose value at distance 3*l

# Then the same but domain size just sligtly larger than 6*l with data in the center

# Then the same but with data at the border

using DIVAnd
using Compat: @info, range
using PyPlot
if VERSION >= v"0.7"
    using LinearAlgebra
end

include("./prep_dirs.jl")

# observations
x = [150.];
y = [150.];
f = [1.];

# final grid
xi,yi = ndgrid(range(1,stop=299,length=299),range(1,stop=299,length=299));

# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);

# correlation length
len = 10;

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field
fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);

@show fi[150,150]
@show fi[150,179]

# observations
x = [30.];
y = [30.];
f = [1.];

# final grid
xi,yi = ndgrid(range(1,stop=59,length=59),range(1,stop=59,length=59));

# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);
#Test to push boundary to wider distance:

alen=1.5*len

pn[:,59]=pn[:,59]/alen;
pn[:,1]=pm[:,1]/alen;
pm[59,:]=pm[59,:]/alen;
pm[1,:]=pm[1,:]/alen;


# correlation length
len = 10;

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field
fi2,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);


pcolor(reshape(diag(s.P),59,59)')
colorbar()

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_1.png")));
savefig(figname)
@info "Saved figure as " * figname

@show fi2[30,30]
@show fi2[30,59]

x = [30.];
y = [58.95];
f = [1.];

# final grid
xi,yi = ndgrid(range(1,stop=59,length=59),range(1,stop=59,length=59));



# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);
alen=1.5*len
pn[:,59]=pn[:,59]/alen;
pn[:,1]=pm[:,1]/alen;
pm[59,:]=pm[59,:]/alen;
pm[1,:]=pm[1,:]/alen;


# correlation length
len = 10;

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field
fi3,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);


@show fi3[30,30]
@show fi3[30,59]

#pcolor(reshape(diag(s.P),59,59))



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
