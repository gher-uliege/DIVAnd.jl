# A simple example of DIVAnd in 1 dimension
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using Compat: @info, range
using PyPlot

include("./prep_dirs.jl")

# observations with some points outside
x = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1]
f = sin.(x*6);

# final grid
xi=collect(range(0,stop=1,length=30));

# reference field
fref = sin.(xi*6) ;

# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
pm = ones(size(xi)) / (xi[2]-xi[1]);

# correlation length
len = 0.1;

# obs. error variance normalized by the background error variance
epsilon2 = 0.1;

# fi is the interpolated field
fi,s = DIVAndrun(mask,(pm,),(xi,),(x,),f,len,epsilon2;alphabc=0);

plot(x,f,".",label="observations")
plot(xi,fi,"-",label="analysis")
legend()
figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_results.png")))
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
