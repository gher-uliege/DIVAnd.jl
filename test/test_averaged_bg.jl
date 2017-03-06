# A simple example of divand in 2 dimensions
# with observations from an analytical function.
using Base.Test

using divand
srand(1)
x=randn(100)
y=randn(100)
z=randn(100)
t=randn(100)
f=z

mask,(pm,pn,po,pq),(xi,yi,zi,ti) = divand_rectdom(linspace(-1,1,5),linspace(-1,1,5),linspace(-1,1,5),linspace(-1,1,5))

# correlation length
len = 1

# obs. error variance normalized by the background error variance
epsilon2 = 0.01;

fi,fanom= divand_averaged_bg(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2,[true true false true]);

@test -1.6 < fi[1,1,1,1] < -1.4


# Copyright (C) 2014, 2017 Alexander Barth <a.barth@ulg.ac.be>
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
