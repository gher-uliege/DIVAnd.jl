# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot
if VERSION >= v"0.7"
    using Random
    using Statistics
    using LinearAlgebra
end


x=[0]
y=[0]
z=[0]
f=[1]

if VERSION >= v"0.7"
   Random.seed!(1)
else
   srand(1)
end
nobs=2
x=randn(nobs)
y=randn(nobs)
z=randn(nobs)
f=x+y.*z

jsize=50
mask,(pm,pn,po),(xi,yi,zi) = DIVAnd_rectdom(range(-1,stop=1,length=jsize),range(-1,stop=1,length=jsize),range(-1,stop=1,length=jsize))

# correlation length
len = (0.4,0.4,0.4)

# obs. error variance normalized by the background error variance
epsilon2 = 0.010;

#@time fi,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2;alphabc=1);

epsilon2b=(1+epsilon2)^0.2-1
#epsilon2b=epsilon2

alpha1D=[]

#epsilon2b=epsilon2
@time fi1,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,(0.4*1,0.,0.),epsilon2b;alphabc=1,alpha=alpha1D);

PC1=s.P
H1=s.H
xg1=statevector_pack(s.sv,(fi1,))




@time fi2,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,(0., 0.4*1,0.),epsilon2b;alphabc=1,alpha=alpha1D);

PC2=s.P
H2=s.H

xg2=statevector_pack(s.sv,(fi2,))


@time fi3,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,(0.,0.,0.4*1),epsilon2b;alphabc=1,alpha=alpha1D);

@show extrema(s.H-H1)
@show extrema(s.H-H2)
@show extrema(H2-H1)



PC3=s.P
xg3=statevector_pack(s.sv,(fi3,))


# try the real R ?
zzz=(s.H'* (s.R \ s.yo[:]))

#(s.H'* (s.R \(s.H*(

@show size(zzz)

xgs=PC1*(s.H'* (s.R \(s.H*(PC2*(s.H'* (s.R \(s.H*(PC3*(s.H'* (s.R \(s.H*(PC2*(s.H'* (s.R \(s.H*(PC1*zzz))))))))))))))))


xgs=PC1*(s.H'* (s.R \(s.H*(PC2*(s.H'* (s.R \(s.H*(PC3*(s.H'* (s.R \(s.H*(PC3*(s.H'* (s.R \(s.H*(PC2*(s.H'* (s.R \(s.H*(PC1*zzz))))))))))))))))))))

xguess=statevector_unpack(s.sv,xgs)
xgs2=PC1*(PC2*(PC3*(PC3*(PC2*(PC1*zzz)))))*epsilon2b^5

xgs3=statevector_pack(s.sv,((fi1+fi2+fi3)/3,))

scalef2=dot(xgs,xgs2)/dot(xgs2,xgs2)
@show scalef2
@show epsilon2^4
scalef2=0.0005

vv=var(scalef2*xgs2-xgs)/var(xgs)
@show vv

#   @show var((fi1+fi2+fi3)/3-fi)/var(fi)

tol = 1e-4


maxiter=1000

pcargs = [(:tol, tol),(:maxit,maxiter)]



diagshift=0.00004;

#PC=(PC1+PC2+PC3)/3.
#PC=PC*PC
xr=randn(jsize^3,1)
scalef=xr'*(PC1*(PC2*(PC3*(PC3*(PC2*(PC1*xr))))))./(xr'*xr)
@show scalef
scalef2=1/scalef[1]
function compPC(iB,H,R)
    function fun!(x,fx)
        fx[:] = diagshift*x+scalef2*(PC1*(PC2*(PC3*(PC3*(PC2*(PC1*x))))))
    end
    return fun!
end
# First guess is the HI* coarse solution

# HI*(sc.P*(HI'  *x ))  should be a good operator for the M-1 x operation in preconditionner ?
# Why do I need to take sc.P\ ??? So better use components of P*HI' ?

sv = s.sv
s = 0

@time fiiter,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2;alphabc=1,pcargs...,inversion=:pcg,compPC = compPC, fi0 = unpack(sv,xgs)[1],btrunc=1);

#@time fiiter2,s = DIVAndrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2;alphabc=1,pcargs...,inversion=:pcg,btrunc=1)#, fi0 =xguess);
# Then run with normal resolution and preconditionner

#var(fi-fiiter)/var(fi)


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
