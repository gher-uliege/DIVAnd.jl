# A simple example of DIVAnd in 1 dimensions
# with observations from an analytical function.

using DIVAnd
using Compat: @info, range
using PyPlot
using Interpolations
if VERSION >= v"0.7"
    using Printf
end

include("./prep_dirs.jl")

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

function detrend(mask,pm,xi,x,f,len,epsilon2;
                 niter = 10,
                 progressiter = (i,fi) -> nothing
                 )

    nlevels = length(mask)

    fi = [zeros(size(xi[i][1])) for i = 1:nlevels]



    tmp = copy(f)
    # iterative solver
    for i = 1:niter

        # loop over all levels starting with the finest
        for k = nlevels:-1:1
            # remove all scales from observations f expect scale k
            tmp[:] = f
            for j = 1:nlevels
                if j != k
                    tmp = f - DIVAnd.interp(xi[j],fi[j],x)
                end
            end

            # analyse variations at scale k
            fi[k][:],s = DIVAndrun(mask[k],pm[k],xi[k],x,tmp,len[k],epsilon2);
        end

        progressiter(i,fi)
    end

    # final analysis:
    # add all scales

    fa = copy(fi[nlevels])

    for i = 1:nlevels-1
        fa = fa + DIVAnd.interp(xi[i],fi[i],xi[nlevels])
    end

    return fa,fi
end



# grid
# xi[i][j][k₁,k₂,...] is the coordinate of the point k₁,k₂,...
# along the dimension j for the i-th grid
#
# The last grid is the final grid

xi = (
    (collect(range(0,stop=4*pi,length=400)),),  # coarse grid
    (collect(range(0,stop=4*pi,length=400)),)   # fine grid
)


nlevels = length(xi)

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

len = [1.5, 1.5/10];

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field at different scales
# fi[1] is u
# fi[2] is x


function plotiter(i,fi)
    figure(1)
    for k = 1:nlevels
        subplot(nlevels,1,k)
        plot(xi[k][1],fi[k],"-",label="iteration $(i)");
        title("level $(k)")
    end

    figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => @sprintf("_%04d.png",i))));
    savefig(figname)
    @info "Saved figure as " * figname

end

fa,fi = detrend(mask,pm,xi,x,f,len,epsilon2;
                niter = 10,
                progressiter = plotiter,
                )

figure(2)
plot(x[1],f,".",label="observation")
plot(xi[1][1],fi[1],"-",label="analysis (trend)")
#plot(xi[2][1],fi[2],"-",label="analysis (variations)")
plot(xi[2][1],fa,"-",label="analysis (total)")
legend()

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_2.png")));
savefig(figname)
@info "Saved figure as " * figname

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
