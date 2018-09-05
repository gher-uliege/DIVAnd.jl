#SBATCH --mem-per-cpu=12000


# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot
if VERSION >= v"0.7"
    using LinearAlgebra
    using Statistics
end


x=[0]
y=[0]
z=[0]
f=[1]

x=randn(200)
y=randn(200)
z=randn(200)
t=randn(200)
f=x+y-t+z

mask,(pm,pn,po,pq),(xi,yi,zi,ti) = DIVAnd_rectdom(range(-1,stop=1,length=15),range(-1,stop=1,length=15),range(-1,stop=1,length=15),range(-1,stop=1,length=15))

# correlation length
len = (0.4,0.4,0.4,0.4)

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

@time fi,s = DIVAndrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2;alphabc=2);

@time fipca,spc = DIVAndrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,(0.4,0.4,0.,0.),epsilon2;alphabc=2);
PCA=spc.P
#mpca=mean(diag(PCA))
#@show mpca
xguessa=statevector_pack(spc.sv,(fipca,))
mga=mean(xguessa)
@show mga
@time fipcb,spcb = DIVAndrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,(0.,0.,0.4,0.4),epsilon2;alphabc=2);
PCB=spcb.P
#mpcb=mean(diag(PCB))
#@show mpcb
xguessb=statevector_pack(spcb.sv,(fipcb,));
mgb=mean(xguessb)
@show mgb
xguess=(xguessa.+xguessb)/2;
tol = 2e-3


maxiter=10000

pcargs = [(:tol, tol),(:maxit,maxiter)]



diagshift=0.004;



function compPC(iB,H,R)
    function fun!(x,fx)
        fx[:] = diagshift*x+((PCA*x)+(PCB*x))/2;
    end
    return fun!
    #return x -> diagshift*x+0.30698675.*(PCA*(PCB*x));
    #return x -> diagshift*x+(PCA*(x-PCB*x)+PCB*(x-PCA*x));
    #return x -> diagshift*x+(PCA*x);
    #     return jmPHI'*(jmPHI*x);
    #   return x->x;
end
# First guess is the HI* coarse solution

# HI*(sc.P*(HI'  *x ))  should be a good operator for the M-1 x operation in preconditionner ?
# Why do I need to take sc.P\ ??? So better use components of P*HI' ?


# Then run with normal resolution and preconditionner
@time fiiter,s = DIVAndrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2;alphabc=2,pcargs...,
                           inversion=:pcg,compPC = compPC, fi0 = unpack(spc.sv,xguess)[1]);

@time fiiter2,s = DIVAndrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2;alphabc=2,pcargs...,
                            inversion=:pcg, fi0 = unpack(spc.sv,xguess)[1]);

var(fi-fiiter)/var(fi)

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
