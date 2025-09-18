# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using LinearAlgebra
using Statistics
using Makie, CairoMakie

include("./prep_dirs.jl")
figbasename = basename(@__FILE__)

# observations
x = [-0.7, 0.7, 0.7, -0.7];
y = [0.0, -0.7, 0.7, 0.7];
f = [1.0, -1.0, 1.0, -1.0];

# final grid
xi, yi = ndgrid(range(-1, stop = 1, length = 81), range(-1, stop = 1, length = 81));
xifin, yifin =
    ndgrid(range(-10, stop = 10, length = 801), range(-10, stop = 10, length = 801));

varb1 = 0
Bi = 0
Bold = 0

pmfin = ones(size(xifin)) / (xifin[2, 1] - xifin[1, 1]);
pnfin = ones(size(xifin)) / (yifin[1, 2] - yifin[1, 1]);

# correlation length
len = 1 / 4.1 * 8 / 9.75;
@show len
@show 2 / len

# obs. error variance normalized by the background error variance
epsilon2 = 1.0;

maskfin = trues(size(xifin));
# fi is the interpolated field
fifin, sfin = DIVAndrun(
    maskfin,
    (pmfin, pnfin),
    (xifin, yifin),
    (x, y),
    f,
    len,
    epsilon2;
    alphabc = 0,
);

firef = fifin[401-40:401+40, 401-40:401+40];


# all points are valid points
mask = trues(size(xi));

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2, 1] - xi[1, 1]);
pn = ones(size(xi)) / (yi[1, 2] - yi[1, 1]);

@show pm[1, 1] * len
# correlation length


# obs. error variance normalized by the background error variance
epsilon2 = 1.0;

# fi is the interpolated field
fi, s = DIVAndrun(mask, (pm, pn), (xi, yi), (x, y), f, len, epsilon2; alphabc = 1.09);

Bnew = diag(s.P)
@show mean(Bnew)

fiold, sold = DIVAndrun(mask, (pm, pn), (xi, yi), (x, y), f, len, 1E6; alphabc = 0);
Bold = diag(sold.P);
@show mean(Bold)

#pcolor(reshape(Bold,81,81));colorbar()

fiold, s = DIVAndrun(mask, (pm, pn), (xi, yi), (x, y), f, len, epsilon2; alphabc = 0);

@show sqrt(var(fiold - fi)) / sqrt(var(fiold - firef))
varr = zeros(100)
rms = zeros(100)
al = zeros(100)
for ii = 1:100
    alen = 0.25 + ii / 100 * 1.25
    fi, si = DIVAndrun(mask, (pm, pn), (xi, yi), (x, y), f, len, epsilon2; alphabc = alen)
    fbi, sbi = DIVAndrun(mask, (pm, pn), (xi, yi), (x, y), f, len, 1E6; alphabc = alen)
    Bi = diag(sbi.P)
    varr[ii] = sqrt(var(Bi)) / sqrt(var(Bold))
    rms[ii] = sqrt(var(firef - fi)) / sqrt(var(fiold - firef))
    al[ii] = alen

end

fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, al, varr)
ax2 = Axis(fig[1, 2])
lines!(ax2, al, rms)

@show al[argmin(varr)]

figname = joinpath(figdir, replace(figbasename, r".jl$" => "_1.png"));
save(figname, fig)
@info "Saved figure as " * figname


alphaopt = al[argmin(varr)]
fi, swhat =
    DIVAndrun(mask, (pm, pn), (xi, yi), (x, y), f, len, epsilon2; alphabc = alphaopt);
#fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;alphabc=0);
fidir, swhat = DIVAndrun(mask, (pm, pn), (xi, yi), (x, y), f, len, epsilon2; alphabc = 2.06);


fig = Figure();

rmval = round(sqrt(var(firef - firef)), digits = 4)
ax1 = Axis(fig[1, 1], title = "Reference rms $rmval")
heatmap!(ax1, xi[:, 1], yi[1, :], firef);
Colorbar(fig[1, 2], limits = (-1, 1))
scatter!(ax1, x, y, color = :black)

rmval = round(sqrt(var(fi - firef)), digits = 4)
ax2 = Axis(fig[1, 3], title = "optimal alpha  rms $rmval")
heatmap!(ax2, xi[:, 1], yi[1, :], fi);
Colorbar(fig[1, 4], limits = (-1, 1))

rmval = round(sqrt(var(fidir - firef)), digits = 4)
ax3 = Axis(fig[2, 1], title = "alpha=1.09  rms $rmval")
heatmap!(ax3, xi[:, 1], yi[1, :], fidir);
Colorbar(fig[2, 2], limits = (-1, 1))

rmval = round(sqrt(var(fiold - firef)), digits = 4)
ax4 = Axis(fig[2, 3], title = "old bc  rms $rmval")
heatmap!(ax4, xi[:, 1], yi[1, :], fiold);
Colorbar(fig[2, 4], limits = (-1, 1))

figname = joinpath(figdir, replace(figbasename, r".jl$" => "_2.png"));
save(figname, fig)
@info "Saved figure as " * figname


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
