# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

if VERSION >= v"0.7"
    using Test
else
    using Base.Test
    using Compat: range
end

using DIVAnd


# observations
x = [0.5];
y = [0.5];
f = [1.];

# final grid
mask,(pm,pn),(xi,yi) = DIVAnd_rectdom(
    range(0,stop=1,length=23),
    range(0,stop=1,length=22))

# correlation length
len = 0.005;

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field

fiex,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,(0.5*len,1*len),epsilon2);

# to compare to the filtered version
fiexf=DIVAnd_filter3(fiex,NaN,2)

fi,s = DIVAndgo(mask,(pm,pn),(xi,yi),(x,y),f,(0.5*len,1*len),epsilon2);

fifp,s = DIVAndgo(mask,(pm,pn),(xi,yi),(x,y),f,(0.5*len,1*len),epsilon2;moddim=[0,0]);

@test maximum(fi) ≈ maximum(fiexf)

@test maximum(fifp) ≈ maximum(fiexf)


# Copyright (C) 2014, 2017 Alexander Barth         <a.barth@ulg.ac.be>
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
