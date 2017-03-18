# A simple example of divand in 2 dimensions
# with observations from an analytical function.

using divand
using PyPlot


x=[0]
y=[0]
z=[0]
f=[1]

srand(876)
nobs=2
x=randn(nobs)
y=randn(nobs)
z=randn(nobs)
t=randn(nobs)
f=x+y.*z+t.*t.*x

@show var(f)

jsize=40
mask,(pm,pn,po,pq),(xi,yi,zi,ti) = divand_rectdom(linspace(-1,1,jsize),linspace(-1,1,jsize),linspace(-1,1,jsize),linspace(-1,1,jsize))

# correlation length
len = (0.4,0.4,0.4,0.4)

# obs. error variance normalized by the background error variance
epsilon2 = 1/4;



epsilon2b=epsilon2*100
@time fi1,s = divandrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,(0.4,0.4,0,0),epsilon2b;alphabc=1,btrunc=2);

PC1=s.P
xg1=statevector_pack(s.sv,(fi1,))
s=0
gc()

@time fi2,s = divandrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,(0,0,0.4,0.4),epsilon2b;alphabc=1,btrunc=2);


PC2=s.P
xg2=statevector_pack(s.sv,(fi2,))
s=0
gc()

   xguess=(xg1.+xg2)/2.;     
   
#   @show var((fi1+fi2+fi3)/3-fi)/var(fi)
   
        tol = 1e-4


        maxiter=1000

        pcargs = [(:tol, tol),(:maxit,maxiter)]



        diagshift=0.00004;

#PC=(PC1+PC2+PC3)/3.
#PC=PC*PC
xr=randn(jsize^4,1)
scalef=xr'*(PC1*(PC2*(((PC2*(PC1*xr))))))./(xr'*xr)
@show scalef

function compPC(iB,H,R)
#            return x -> diagshift*x+PC*x;
#            return x -> diagshift*x+1./9.*PC1*(PC1*x+PC2*x+PC3*x)+1./9.*PC2*(PC1*x+PC2*x+PC3*x)+1./9.*PC3*(PC1*x+PC2*x+PC3*x)			    
#return x -> diagshift*x+1./9.*(PC1*(PC1*x+(PC2*x+(PC3*x))))+1./9.*(PC2*(PC1*x+(PC2*x+(PC3*x))))+1./9.*(PC3*(PC1*x+(PC2*x+(PC3*x))))   	
#return x -> diagshift*x+1./3.*(PC1*(PC1*x)+(PC2*(PC2*x)+(PC3*(PC3*x))))#+1./9.*(PC2*(PC1*x+(PC2*x+(PC3*x))))+1./9.*(PC3*(PC1*x+(PC2*x+(PC3*x))))   	
      return x -> diagshift*x+1./scalef[1]*(PC1*(PC2*(((PC2*(PC1*x))))))

			#return x -> diagshift*x+0.30698675.*(PCA*(PCB*x));
			#return x -> diagshift*x+(PCA*(x-PCB*x)+PCB*(x-PCA*x));
			#return x -> diagshift*x+(PCA*x);
            #     return jmPHI'*(jmPHI*x);
            #   return x->x;
        end
        # First guess is the HI* coarse solution

        # HI*(sc.P*(HI'  *x ))  should be a good operator for the M-1 x operation in preconditionner ?
# Why do I need to take sc.P\ ??? So better use components of P*HI' ?,ti


@time fiiter,s = divandrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2;alphabc=1,pcargs...,inversion=:pcg,compPC = compPC,btrunc=2)#, fi0 =xguess);

@time fiiter2,s = divandrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2;alphabc=1,pcargs...,inversion=:pcg,btrunc=2)#, fi0 =xguess);
# Then run with normal resolution and preconditionner





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
