#SBATCH --mem-per-cpu=16000


# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot
if VERSION >= v"0.7"
    using Random
    using Statistics
end


x=[0.]
y=[0.]
z=[0.]
f=[1.]

if VERSION >= v"0.7.0-beta.0"
   Random.seed!(876)
else
   srand(876)
end
nobs=100
x=randn(nobs)
y=randn(nobs)
z=randn(nobs)
t=randn(nobs)
f=x+y.*z+t.*t.*x

@show var(f)

jsize=100
jsizez=7
jsizet=12

# jsize=40
# jsizez=40
# jsizet=40
mask,(pm,pn,po,pq),(xi,yi,zi,ti) = DIVAnd_rectdom(range(-1,stop=1,length=jsize),range(-1,stop=1,length=jsize),range(-1,stop=1,length=jsizez),range(-1,stop=1,length=jsizet))

# correlation length
#len = (0.4,0.4,0.2,0.2)

len = (0.3,0.2,0.1,0.1)

# obs. error variance normalized by the background error variance
epsilon2 = 0.3;



epsilon2b=1000.
@time fi1,s = DIVAndrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,(0.,0.,len[3]/1.42,len[4]/1.42),epsilon2b;alphabc=1);

sv = s.sv
PC1=s.P
xg1=statevector_pack(s.sv,(fi1,))
s=0

@show size(PC1)

@time fi2,s = DIVAndrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,(len[1],len[1],0.,0.),epsilon2;alphabc=1);


PC2=s.P
xg2=statevector_pack(s.sv,(fi2,))
s=0

xguess=PC1*(PC1*xg2);
xr=randn(jsize^2*jsizez*jsizet,1)
scaleff=xr'*(PC1*(PC1*xr))./(xr'*xr)
@show scaleff


#   @show var((fi1+fi2+fi3)/3-fi)/var(fi)

tol = 1e-4


maxiter=10*Int(ceil(sqrt(size(xg2)[1])))
@show maxiter

pcargs = [(:tol, tol),(:maxit,maxiter)]


diagshift=0.00004;
diagshift=0.00001;

#PC=(PC1+PC2+PC3)/3.
#PC=PC*PC

scalef=xr'*(PC1*(PC2*((((PC1*xr))))))./(xr'*xr)

scalefbis=xr'*(PC1*(((((PC1*xr))))))./(xr'*xr)
scalefter=xr'*((((((PC2*xr))))))./(xr'*xr)
@show scalef,scalefbis,scalefter

@show scalef
scalef2=1/scalef[1]
@show scalef2

# This strange estimation worked quite well
scalef2=2*scalef2/sqrt(nobs)*sqrt(epsilon2)

@show scalef2,1/scaleff[1]
# Other test ...
scalef2=scalefter[1]/scalef[1]
xguess=xguess*scalef2
function compPC(iB,H,R)
    #            return x -> diagshift*x+PC*x;
    #            return x -> diagshift*x+1./9.*PC1*(PC1*x+PC2*x+PC3*x)+1./9.*PC2*(PC1*x+PC2*x+PC3*x)+1./9.*PC3*(PC1*x+PC2*x+PC3*x)
    #return x -> diagshift*x+1./9.*(PC1*(PC1*x+(PC2*x+(PC3*x))))+1./9.*(PC2*(PC1*x+(PC2*x+(PC3*x))))+1./9.*(PC3*(PC1*x+(PC2*x+(PC3*x))))
    #return x -> diagshift*x+1./3.*(PC1*(PC1*x)+(PC2*(PC2*x)+(PC3*(PC3*x))))#+1./9.*(PC2*(PC1*x+(PC2*x+(PC3*x))))+1./9.*(PC3*(PC1*x+(PC2*x+(PC3*x))))
    function fun!(x,fx)
        fx[:] = diagshift*x+scalef2*(PC1*(PC2*((((PC1*x))))))
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
# Why do I need to take sc.P\ ??? So better use components of P*HI' ?,ti

@time fiiter,s = DIVAndrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2;alphabc=1,pcargs...,
                           inversion=:pcg,compPC = compPC,btrunc=2, fi0 = unpack(sv,xguess)[1]);

#@time fiiter2,s = DIVAndrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2;alphabc=1,pcargs...,inversion=:pcg,btrunc=2)#, fi0 =xguess);
# Then run with normal resolution and preconditionner



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
