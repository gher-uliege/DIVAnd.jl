# A simple example of divand in 1 dimensions
# with observations from an analytical function.

using divand
using PyPlot
using Interpolations

function interp!(xi::NTuple{N,Array{T,N}},fi::Array{T,N},x,f) where T where N
    # tuple of vector with the varying parts
    xivector = ntuple(j -> xi[j][[(i==j ? (:) : 1 ) for i in 1:N]...], N) :: NTuple{N,Vector{T}}
    
    itp = interpolate(xivector,fi,Gridded(Linear()))

    xpos = zeros(N)
    for i in eachindex(f)
        # position of the i-th location in f
        for j = 1:N
            xpos[j] = x[j][i]
        end
        f[i] = itp[xpos...]
    end
end

"""
    f = interp(xi,fi,x)

Interpolate field `fi` (n-dimensional array) defined at `xi` (tuble of 
n-dimensional arrays) onto grid `x` (tuble of n-dimensional arrays).
The grid in `xi` must be align with the axis (e.g. produced by ndgrid).
"""

function interp(xi,fi,x)
    # check size
    @assert all([size(xc) == size(fi) for xc in xi])
    
    f = similar(x[1])
    interp!(xi,fi,x,f)    
    return f
end

    

# grid
# xi[i][j][k₁,k₂,...] is the coordinate of the point k₁,k₂,...
# along the dimension j for the i-th grid
#
# The last grid is the final grid

xi = (
    (collect(linspace(0,4*pi,200)),),  # coarse grid
    (collect(linspace(0,4*pi,200)),)   # fine grid
)


nlevels = 2

# sample data on fine grid
fitrue = sin.(xi[2][1]) + sin.(10*xi[2][1])


# observations
ind = 1:2:length(xi[2][1])
x = (xi[2][1][ind],)
f = fitrue[ind];

# all points are valid points
mask = ntuple(i -> trues(size(xi[i][1])), length(xi));

# this problem has a simple cartesian metric
# pm[1] is the inverse of the resolution along the 1st dimension of grid 1
# pm[2] is the inverse of the resolution along the 1st dimension of grid 2

pm = ntuple(i -> (ones(size(xi[i][1])) / (xi[i][1][2]-xi[i][1][1]),), length(xi))


# correlation length for different scales
# len[1] -> long-term trend (u)
# len[2] -> short-term variation (x)

len = [1, 0.1];

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field at different scales
# fi[1] is u
# fi[2] is x

fi = [zeros(xi[i][1]) for i = 1:2]

niter = 100

# iterative solver
for i = 1:niter
    # Anomalies d̃ = d - H₂ u used in classical analysis J₁ and provides x
    
    # remove large-scale trend from f
    tmp = f - interp(xi[1],fi[1],x)
    # find short-term variation x
    fi[2][:],s = divandrun(mask[2],pm[2],xi[2],x,tmp,len[2],epsilon2);

    
    # Anomalies d̃ = d - H₁ x used in classical analysis J₂ and provides u

    # remove short-term variation x from f
    tmp = f - interp(xi[2],fi[2],x)
    # find large-scale trend u
    fi[1][:],s = divandrun(mask[1],pm[1],xi[1],x,tmp,len[1],epsilon2);

    figure(1)
    subplot(2,1,1)
    plot(xi[1][1],fi[1],"-",label="analysis (trend) $(i)");
    title("trend")
    
    subplot(2,1,2)    
    plot(xi[2][1],fi[2],"-",label="analysis (variations) $(i)");
    title("variations")
    
end

# final analysis
# add all parts

fa = copy(fi[nlevels])

for i = 1:nlevels-1
    fa = fa + interp(xi[i],fi[i],xi[nlevels])
end




figure(2)
plot(x[1],f,".",label="observation")
plot(xi[1][1],fi[1],"-",label="analysis (trend)")
#plot(xi[2][1],fi[2],"-",label="analysis (variations)")
plot(xi[2][1],fa,"-",label="analysis (total)")
legend()



# Copyright (C) 2018 Jean-Marie Beckers <jm.beckers@ulg.ac.be>
#               2018 Alexander Barth <a.barth@ulg.ac.be>
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
