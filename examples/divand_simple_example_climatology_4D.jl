# A simple example of divand in 4 dimensions
# with observations from an analytical function.

using divand
using PyPlot
include("../src/override_ssmult.jl")
# final grid
#
testsizex=100
testsizey=90
testsizez=10
testsizet=12
# observations
nobs=2000;
x = rand(nobs)*testsizex;
y = rand(nobs)*testsizey;
z = rand(nobs)*testsizez;
t = rand(nobs)*testsizet;
f = sin.(x*pi/180) .* cos.(y*pi/180.)+sin.(z*6/50) .* cos.(x*6*pi/180) .* sin.(t*2*pi/12);


xi,yi,zi,ti = ndgrid(linspace(1,testsizex,testsizex),linspace(1,testsizey,testsizey),linspace(1,testsizez,testsizez),linspace(1,testsizet,testsizet));

# reference field
fref = sin.(xi*pi/180) .* cos.(yi*pi/180.)+sin.(zi*6/50) .* cos.(xi*6*pi/180) .* sin.(ti*2*pi/12);

# all points are valid points
mask = trues(xi);

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(xi) / (xi[2,1,1,1]-xi[1,1,1,1]);
pn = ones(xi) / (yi[1,2,1,1]-yi[1,1,1,1]);
po = ones(xi) / (zi[1,1,2,1]-zi[1,1,1,1]);
pq = ones(xi) / (ti[1,1,1,2]-ti[1,1,1,1]);

# correlation length
len = (8, 8, 0, 0);

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field
@time fi,s = divandrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2; moddim=[0,0,0,12]);

PP=s.P
x0=statevector_pack(s.sv,(fi,))

s=0
fi=0

gc()

len = (8, 8, 1, 1);

# obs. error variance normalized by the background error variance


tol = 5e-2


maxiter=100
@show maxiter

pcargs = [(:tol, tol),(:maxit,maxiter)]



diagshift=0.001;
@show diagshift

function compPC(iB,H,R)
    return x -> diagshift*x+PP*x;
    #     return jmPHI'*(jmPHI*x);
    #   return x->x;
end


# fi is the interpolated field
@time fi,s = divandrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2; moddim=[0,0,0,12], pcargs..., inversion=:pcg,compPC = compPC, fi0 =x0);


# Copyright (C) 2014, 2017 Alexander Barth <a.barth@ulg.ac.be>
#                         Jean-Marie Beckers   <JM.Beckers@ulg.ac.be>
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
