# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot


x=[0.]
y=[0.]
z=[0.]
f=[1.]

x=randn(200)
y=randn(200)
z=randn(200)
f=x+y+z

mask,(pm,pn,po),(xi,yi,zi) = DIVAnd_rectdom(range(-1,stop=1,length=30),range(-1,stop=1,length=30),range(-1,stop=1,length=30))

# correlation length
len = 0.3

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

@time fipca,spc = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2;alphabc=1,alpha=[1,2,1]);
PCA=spc.P
iB2=spc.iB

xguess=fipca


tol = 2e-3


maxiter=10000

pcargs = [(:tol, tol),(:maxit,maxiter)]



diagshift=0.0004;



function compPC(iB,H,R)
    function fun!(x,fx)
        fx[:] = diagshift*x+  PCA*x .- 0.5.*(PCA*(iB*(PCA*x)-iB2*(PCA*x)));
    end
    return fun!
end




@time fiiter,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2;alphabc=1,pcargs...,inversion=:pcg,compPC = compPC, fi0 =xguess);

nothing

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
