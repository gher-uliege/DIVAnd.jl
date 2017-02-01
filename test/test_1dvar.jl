# Testing divand in 1 dimension.

# grid of background field
xi = collect(linspace(0,1,21));

x = [.4; .6];
f = [.4; .6];

mask = trues(xi);
mask[[1 end]] = false;
  
pm = ones(xi) / (xi[2]-xi[1]);

len = 0.1
epsilon2 = 0.5
  
fi,s = divandrun(mask,(pm,),(xi,),(x,),f,len,epsilon2);

@test xi[fi .== maximum(fi)][1] == x[2]


# Copyright (C) 2014,2017 Alexander Barth <a.barth@ulg.ac.be>
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
