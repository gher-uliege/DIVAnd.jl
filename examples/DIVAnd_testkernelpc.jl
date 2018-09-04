# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot

if VERSION >= v"0.7"
    using Statistics
end

x=[0.]
y=[0.]
z=[0.]
f=[1.]

x=randn(200)
y=randn(200)
z=randn(200)
f=x+y

mask,(pm,pn,po),(xi,yi,zi) = DIVAnd_rectdom(range(-1,stop=1,length=40),range(-1,stop=1,length=40),range(-1,stop=1,length=40))

# correlation length
len = (0.2,0.2,0.2)

# obs. error variance normalized by the background error variance
epsilon2 = 2.;

@time fi,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2;alphabc=2);

itest=20
jtest=20
rms=zeros(itest)
epsfac=zeros(itest)
rms2=zeros(itest)
epsfac2=zeros(itest)

ffac=zeros(itest)
ffac2=zeros(itest)
rmsb=zeros(itest)
rmsb2=zeros(itest)


lenf=zeros(itest)
fipc=0
fipc2=0

ii=1

lenfac=0.5+ii/itest*4.0

lenfac=4.1
epsfacc=0.00005+(ii-1)/itest*0.02



lena = ([x[1]*lenfac for x in len]...,)

@time fipc2,spc2 = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,lena,epsilon2*epsfacc;alphabc=2,alpha=[1,1]);

xguess=fipc2;
scP=spc2.P;


tol = 2e-3


maxiter=10000

pcargs = [(:tol, tol),(:maxit,maxiter)]



diagshift=0.004;



function compPC(iB,H,R)
    function fun!(x,fx)
        fx[:] = diagshift*x+(scP*x)
    end
    return fun!
end
# First guess is the HI* coarse solution

# HI*(sc.P*(HI'  *x ))  should be a good operator for the M-1 x operation in preconditionner ?
# Why do I need to take sc.P\ ??? So better use components of P*HI' ?


# Then run with normal resolution and preconditionner
@time fiiter,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2;alphabc=2,pcargs...,inversion=:pcg,compPC = compPC, fi0 =xguess);

@time fiiter2,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2;alphabc=2,pcargs...,inversion=:pcg, fi0 =xguess);

var(fi-fiiter)/var(fi)

var(fi-fipc2)/var(fi)

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
