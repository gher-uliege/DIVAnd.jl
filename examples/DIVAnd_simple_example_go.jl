#SBATCH --mem-per-cpu=16000


# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

using DIVAnd
using Makie, CairoMakie
using LinearAlgebra
using Random
using Statistics

include("./prep_dirs.jl")
figbasename = basename(@__FILE__)

Random.seed!(1234)

# observations
nobs = 100
x = 0.01 .+ 0.98 * rand(nobs);
y = 0.01 .+ 0.98 * rand(nobs);
x = -0.1 .+ 1.2 * rand(nobs);
y = -0.1 .+ 1.2 * rand(nobs);
f = sin.(x * 6) .* cos.(y * 6);
#f=-1+2*x
# x=[0.5,0.75]
# y=[0.5,0.75]
# f=[1,1]

# final grid
xi, yi = ndgrid(range(0, stop = 1, length = 950), range(0, stop = 1, length = 830));

# reference field
fref = sin.(xi * 6) .* cos.(yi * 6);

# all points are valid points
mask = trues(size(xi));

mask[300:600, 400:600] .= false

# this problem has a simple cartesian metric
# pm is the inverse of the resolution along the 1st dimension
# pn is the inverse of the resolution along the 2nd dimension

pm = ones(size(xi)) / (xi[2, 1] - xi[1, 1]);
pn = ones(size(xi)) / (yi[1, 2] - yi[1, 1]);

# correlation length
len = 0.03;

# obs. error variance normalized by the background error variance
epsilon2 = 1.0;
vscale = 0.001
vscale = 0

# fi is the interpolated field
@time fiexOLD, s =
    DIVAndrun(mask, (pm, pn), (xi, yi), (x, y), f, (len, 0.5 * len), epsilon2; alphabc = 0);

residue = DIVAnd_residualobs(s, fiexOLD)

@time fiOLD, errOLD, residueGO =
    DIVAndgo(mask, (pm, pn), (xi, yi), (x, y), f, (len, 0.5 * len), epsilon2; alphabc = 0);

@show size(residue), size(residueGO)
@show norm(residue - residueGO)

fig = Figure(size = (900, 300))
ax1 = Axis(fig[1, 1], title = "Old version")
heatmap!(ax1, xi[:, 1], yi[1, :], fiOLD)
Colorbar(fig[1, 2], limits = (-1, 1))

ax2 = Axis(fig[1, 3], title = "Reference")
heatmap!(ax2, xi[:, 1], yi[1, :], fiexOLD)
Colorbar(fig[1, 4], limits = (-1, 1))

rms = sqrt(var(fiexOLD - fiOLD))
ax3 = Axis(fig[1, 5], title = "rms $rms")
heatmap!(ax3, xi[:, 1], yi[1, :], fiexOLD - fiOLD)
Colorbar(fig[1, 6], limits = (-0.05, 0.05))
scatter!(ax3, x, y, color = :black)

figname = joinpath(figdir, replace(figbasename, r".jl$" => "_1.png"));
save(figname, fig)
@info "Saved figure as " * figname

@time fiex, s = DIVAndrun(
    mask,
    (pm, pn),
    (xi, yi),
    (x, y),
    f,
    (len, 0.5 * len),
    epsilon2;
    alphabc = 1.0,
);

residue = DIVAnd_residualobs(s, fiex)

@time fi, erri, residueGO =
    DIVAndgo(mask, (pm, pn), (xi, yi), (x, y), f, (len, 0.5 * len), epsilon2; alphabc = 1.0);

@show size(residue), size(residueGO)

@show norm(residue - residueGO)

fig = Figure(size = (900, 300))
ax1 = Axis(fig[1, 1], title = "NEW version")
heatmap!(ax1, xi[:, 1], yi[1, :], fi, colorrange = (-1, 1))

ax2 = Axis(fig[1, 2], title = "Reference")
heatmap!(ax2, xi[:, 1], yi[1, :], fiex, colorrange = (-1, 1))

rms = sqrt(var(fiex - fi))
ax3 = Axis(fig[1, 3], title = "rms $rms")
heatmap!(ax3, xi[:, 1], yi[1, :], fiex - fi, colorrange = (-0.05, 0.05))
scatter!(ax3, x, y, color = :black);
Colorbar(fig[1, 4])

figname = joinpath(figdir, replace(figbasename, r".jl$" => "_2.png"));
save(figname, fig)
@info "Saved figure as " * figname

fig = Figure(size = (800, 400))
ax1 = Axis(fig[1, 1], title = "cpme NEW version")
heatmap!(ax1, xi[:, 1], yi[1, :], erri, colorrange = (0, 1))

ax2 = Axis(fig[1, 2], title = "cpme OLD version")
heatmap!(ax2, xi[:, 1], yi[1, :], errOLD, colorrange = (0, 1))
Colorbar(fig[1, 3])
figname = joinpath(figdir, replace(figbasename, r".jl$" => "_3.png"));
save(figname, fig)
@info "Saved figure as " * figname

fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, residue, residueGO)

# Copyright (C) 2014, 2018 Alexander Barth         <a.barth@ulg.ac.be>
#                          Jean-Marie Beckers   <JM.Beckers@ulg.ac.be>
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
