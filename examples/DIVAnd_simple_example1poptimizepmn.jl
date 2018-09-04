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
    using Statistics
end


len = 4;
aj=zeros(1000)
vj=zeros(1000)

for j=1:1000
    alen=j/100

    # observations
    x = [10.];
    y = [10.];
    f = [1.];

    idim=59
    # final grid
    xi,yi = ndgrid(range(0,stop=30,length=idim),range(0,stop=30,length=idim));



    # all points are valid points
    mask = trues(size(xi));

    # this problem has a simple cartesian metric
    # pm is the inverse of the resolution along the 1st dimension
    # pn is the inverse of the resolution along the 2nd dimension

    pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
    pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);
    #Test to push boundary to wider distance:

    pn[:,idim] .= 1. / (alen*len);
    pn[:,1] .= 1. / (alen*len);
    pm[idim,:] .= 1. / (alen*len);
    pm[1,:] .= 1. / (alen*len);


    # correlation length


    # obs. error variance normalized by the background error variance
    epsilon2 = 10000.;

    # fi is the interpolated field
    fi2,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);


    #pcolor(reshape(diag(s.P),59,59)')
    #colorbar()

    aj[j]=alen
    vj[j]=var(diag(s.P))

end


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
