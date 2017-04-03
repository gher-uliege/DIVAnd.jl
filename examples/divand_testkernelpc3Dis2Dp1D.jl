# A simple example of divand in 2 dimensions
# with observations from an analytical function.

using divand
using PyPlot


x=[0]
y=[0]
z=[0]
f=[1]

srand(1)
nobs=400
x=randn(nobs)
y=randn(nobs)
z=randn(nobs)
f=x+y.*z

jsizex=200*2
jsizey=70*2
jsizet=20

ffull=false
mask,(pm,pn,po),(xi,yi,zi) = divand_rectdom(linspace(-1,1,jsizex),linspace(-1,1,jsizey),linspace(-1,1,jsizet))

# correlation length
len = (0.2,0.2,0.1)

# obs. error variance normalized by the background error variance
epsilon2 = 1;

if ffull
@time fi,s = divandrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2;alphabc=1);
else
fi=zeros(jsizex,jsizey,jsizet)
end

fis=reshape(fi,(jsizex*jsizey*jsizet))

epsilon2b=1000
#epsilon2b=epsilon2

alpha1D=[]

#epsilon2b=epsilon2
@time fi1,s = divandrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,(0,0,len[3]/1.7),epsilon2b;alphabc=1,alpha=alpha1D);

PC1=s.P
H1=s.H
xg1=statevector_pack(s.sv,(fi1,))



@time fi3,s = divandrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,(len[1],len[2],0),epsilon2;alphabc=1,alpha=alpha1D);



PC3=s.P
xg3=statevector_pack(s.sv,(fi3,))





#   @show var((fi1+fi2+fi3)/3-fi)/var(fi)

tol = 1e-4


maxiter=10*Int(ceil(sqrt(jsizex*jsizey*jsizet)))

pcargs = [(:tol, tol),(:maxit,maxiter)]



diagshift=0.00004;
diagshift=0.00001

#PC=(PC1+PC2+PC3)/3.
#PC=PC*PC
xr=randn(jsizex*jsizey*jsizet,1)
scalef=xr'*(PC1*((PC3*(((PC1*xr))))))./(xr'*xr)
@show scalef
scalef2=1/scalef[1]
@show scalef2
scalef2=2*scalef2/sqrt(nobs)*sqrt(epsilon2)



xr=randn(jsizex*jsizey*jsizet,1)
scalef=xr'*(PC1*(PC3*((((PC1*xr))))))./(xr'*xr)

scalefbis=xr'*(PC1*(((((PC1*xr))))))./(xr'*xr)
scalefter=xr'*((((((PC3*xr))))))./(xr'*xr)
@show scalef,scalefbis,scalefter

@show scalef
scalef2=1/scalef[1]
@show scalef2

# This strange estimation worked quite well
scalef2=2*scalef2/sqrt(nobs)*sqrt(epsilon2)


# Other test ...
scalef2=scalefter[1]/scalef[1]
xguess=xg3*scalef2

function compPC(iB,H,R)
    #            return x -> diagshift*x+PC*x;
    #            return x -> diagshift*x+1./9.*PC1*(PC1*x+PC2*x+PC3*x)+1./9.*PC2*(PC1*x+PC2*x+PC3*x)+1./9.*PC3*(PC1*x+PC2*x+PC3*x)
    #return x -> diagshift*x+1./9.*(PC1*(PC1*x+(PC2*x+(PC3*x))))+1./9.*(PC2*(PC1*x+(PC2*x+(PC3*x))))+1./9.*(PC3*(PC1*x+(PC2*x+(PC3*x))))
    #return x -> diagshift*x+1./3.*(PC1*(PC1*x)+(PC2*(PC2*x)+(PC3*(PC3*x))))#+1./9.*(PC2*(PC1*x+(PC2*x+(PC3*x))))+1./9.*(PC3*(PC1*x+(PC2*x+(PC3*x))))
    return x -> diagshift*x+scalef2*(PC1*((PC3*(((PC1*x))))))

    #return x -> diagshift*x+0.30698675.*(PCA*(PCB*x));
    #return x -> diagshift*x+(PCA*(x-PCB*x)+PCB*(x-PCA*x));
    #return x -> diagshift*x+(PCA*x);
    #     return jmPHI'*(jmPHI*x);
    #   return x->x;
end
# First guess is the HI* coarse solution

# HI*(sc.P*(HI'  *x ))  should be a good operator for the M-1 x operation in preconditionner ?
# Why do I need to take sc.P\ ??? So better use components of P*HI' ?

s=0
gc()

@time fiiter,s = divandrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2;alphabc=1,pcargs...,inversion=:pcg,compPC = compPC, fi0 =xguess,btrunc=2);

#@time fiiter2,s = divandrun(mask,(pm,pn,po),(xi,yi,zi),(x,y,z),f,len,epsilon2;alphabc=1,pcargs...,inversion=:pcg,btrunc=1)#, fi0 =xguess);
# Then run with normal resolution and preconditionner

#var(fi-fiiter)/var(fi)



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
