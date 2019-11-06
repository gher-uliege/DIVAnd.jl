# Testing DIVAnd in 3 dimensions without correlation in the 3rd dimension
# (vertically stacked).

using DIVAnd
using Test
using Statistics

# function to interpolate
fun(x, y, z) = sin(6x) * cos(6y) * sin(6z)

# grid of background field

gridsize = (30, 30, 10)
mask, (pm, pn, po), (xi, yi, zi) = DIVAnd_rectdom([range(0, stop = 1, length = s) for s in gridsize]...)

fi_ref = fun.(xi, yi, zi)

ϵ = eps()
# grid of observations
x, y, z = ndgrid(
    range(ϵ, stop = 1 - ϵ, length = 10),
    range(ϵ, stop = 1 - ϵ, length = 10),
    range(ϵ, stop = 1 - ϵ, length = 10),
);

# observations
f = fun.(x, y, z)

# correlation length
len = (fill(0.1, size(mask)), fill(0.1, size(mask)), fill(0., size(mask)))

# obs. error variance normalized by the background error variance
epsilon2 = 0.01;


for alphabc in [0, 1]
    # fi is the interpolated field
    fi, s = DIVAndrun(
        mask,
        (pm, pn, po),
        (xi, yi, zi),
        (x[:], y[:], z[:]),
        f[:],
        len,
        epsilon2;
        alphabc = alphabc,
    )

    fi2 = zeros(size(fi))
    for k = 1:size(mask, 3)
        fi2[:, :, k], s = DIVAndrun(
            mask[:, :, k],
            (pm[:, :, k], pn[:, :, k]),
            (xi[:, :, k], yi[:, :, k]),
            (x[:, :, k][:], y[:, :, k][:]),
            f[:, :, k][:],
            (len[1][:, :, k], len[2][:, :, k]),
            epsilon2;
            alphabc = alphabc,
        )
    end
    @test fi ≈ fi2
end


# Copyright (C) 2014,2017 Alexander Barth <a.barth@ulg.ac.be>
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
