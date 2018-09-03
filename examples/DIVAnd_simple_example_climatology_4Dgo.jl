#SBATCH --mem-per-cpu=62000


# A simple example of DIVAnd in 4 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot
include("../src/override_ssmult.jl")
include("./prep_dirs.jl")

# final grid
#
testsizex=80
testsizey=50
testsizez=5
testsizet=12
# observations
nobs=2000;
x = rand(nobs)*testsizex;
y = rand(nobs)*testsizey;
z = rand(nobs)*testsizez;
t = rand(nobs)*testsizet;
f = sin.(x*pi/180) .* cos.(y*pi/180.)+sin.(z*6/50) .* cos.(x*6*pi/180) .* sin.(t*2*pi/12);


xi,yi,zi,ti = ndgrid(range(1,stop=testsizex,length=testsizex),range(1,stop=testsizey,length=testsizey),range(1,stop=testsizez,length=testsizez),range(1,stop=testsizet,length=testsizet));

# reference field
fref = sin.(xi*pi/180) .* cos.(yi*pi/180.)+sin.(zi*6/50) .* cos.(xi*6*pi/180) .* sin.(ti*2*pi/12);

# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2,1,1,1]-xi[1,1,1,1]);
pn = ones(size(xi)) / (yi[1,2,1,1]-yi[1,1,1,1]);
po = ones(size(xi)) / (zi[1,1,2,1]-zi[1,1,1,1]);
pq = ones(size(xi)) / (ti[1,1,1,2]-ti[1,1,1,1]);

# correlation length
len = (8, 8, 1, 1);

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field

#LPMNR=DIVAnd_Lpmnrange((pm,pn,po,pq),len)
#windowlist,csteps,lmask = DIVAnd_cutter(LPMNR,size(mask),[0,0,0,12])


@time fi,s = DIVAndgo(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2; moddim=[0,0,0,12]);

outputfile = joinpath(outputdir,basename(replace(@__FILE__,r".jl$" => ".nc")));
DIVAnd_save(outputfile,mask,"analysis",fi)
@info "Results written in $outputfile"
# Copyright (C) 2014, 2018 Alexander Barth <a.barth@ulg.ac.be>
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
