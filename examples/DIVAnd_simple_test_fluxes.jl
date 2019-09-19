# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot
if VERSION >= v"0.7"
    using Statistics
end

# function to interpolate
fun(x,y) = sin.(6x) * cos.(6y)

# observations

x = rand(1);
y = rand(1);
f = fun.(x,y)

# final grid
xi,yi = ndgrid(range(0,stop=100,length=100),range(0,stop=110,length=110));

# reference field
fref = fun.(xi,yi)

# all points are valid points
mask = trues(size(xi));

mask[1,:] .= false;
mask[end,:] .= false;
mask[:,1] .= false;
mask[:,end] .= false;

mask[30:80,30:80] .= false;

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);

# correlation length
len = 10;

# obs. error variance normalized by the background error variance
epsilon2 = 10000000.;

h=xi .* (100 .- xi) .+ 20

# fi is the interpolated field

fluxes1=sin.(yi[1,:]./10.)+0.1*rand(size(h)[2])

fluxes2=sin.(xi[:,1]./10.)+0.1*rand(size(h)[1])
rfluxes=0.00000001
rfluxes=1
@time fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;topographyforfluxes=(h,0),fluxes=(fluxes1,0),epsfluxes=rfluxes,alphabc=1,alpha=[1, 0, 1]);

# plotting of results



figure()

pcolor(xi,yi,fi);
colorbar()

title("Interpolated field");

savefig("DIVAnd_simple_example_fluxes1.png")

fluxesafter=zeros(size(h)[2])

for j=1:size(h)[2]
 for i=2:size(h)[1]-2
	if mask[i,j]&& mask[i+1,j]
 		fluxesafter[j]=fluxesafter[j]+h[i,j]*(fi[i+1,j]-fi[i,j])
	end
 end
end
@show var(fluxes1+fluxesafter)
@show var(fluxes1)

@time fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;topographyforfluxes=(0,h),fluxes=(0,fluxes2),epsfluxes=rfluxes,alphabc=1,alpha=[1, 0, 1]);

# plotting of results


figure()

pcolor(xi,yi,fi);
colorbar()

title("Interpolated field");

savefig("DIVAnd_simple_example_fluxes2.png")

fluxesafter=zeros(size(h)[1])

for i=1:size(h)[1]
 for j=2:size(h)[2]-2
	if mask[i,j]&& mask[i,j+1]
 		fluxesafter[i]=fluxesafter[i]+h[i,j]*(fi[i,j+1]-fi[i,j])
	end
 end
end

@show var(fluxes2+fluxesafter)
@show var(fluxes2)


# finally both directions

@time fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;topographyforfluxes=(h,h),fluxes=(fluxes1,fluxes2),epsfluxes=rfluxes,alphabc=1,alpha=[1, 0, 1]);

# plotting of results


figure()
pcolor(xi,yi,fi);
colorbar()

title("Interpolated field");

savefig("DIVAnd_simple_example_fluxes12.png")




fluxesafter=zeros(size(h)[2])

for j=1:size(h)[2]
 for i=2:size(h)[1]-2
	if mask[i,j]&& mask[i+1,j]
 		fluxesafter[j]=fluxesafter[j]+h[i,j]*(fi[i+1,j]-fi[i,j])
	end
 end
end
@show var(fluxes1+fluxesafter)
@show var(fluxes1)

fluxesafter=zeros(size(h)[1])

for i=1:size(h)[1]
 for j=2:size(h)[2]-2
	if mask[i,j]&& mask[i,j+1]
 		fluxesafter[i]=fluxesafter[i]+h[i,j]*(fi[i,j+1]-fi[i,j])
	end
 end
end

@show var(fluxes2+fluxesafter)
@show var(fluxes2)


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
#;
# You should have received a copy of the GNU General Public License along with
# this program; if not, see <http://www.gnu.org/licenses/>.
