# A simple example of divand in 2 dimensions
# with observations from an analytical function.

using divand
using PyPlot


x=[0]
y=[0]
z=[0]
f=[1]

x=randn(200)
y=randn(200)
z=randn(200)
t=randn(200)
f=x+y-t+z

mask,(pm,pn,po,pq),(xi,yi,zi,ti) = divand_rectdom(linspace(-1,1,31),linspace(-1,1,31),linspace(-1,1,31),linspace(-1,1,31))

# correlation length
len = 0.3

# obs. error variance normalized by the background error variance
epsilon2 = 1;

scalel=1.25715/0.69315

#@time fi,s = divandrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2;alphabc=2);

@time fipca,spc = divandjog(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2,[2 2 2 2],[scalel scalel scalel 0];alphabc=2);

       tol = 2e-3


        maxiter=10000

        pcargs = [(:tol, tol),(:maxit,maxiter)]



        diagshift=0.00004;




@time fiiter,s = divandrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2;alphabc=2,pcargs...,inversion=:pcg);

var(fipca-fiiter)/var(fipca)



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
