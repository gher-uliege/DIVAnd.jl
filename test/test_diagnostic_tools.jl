# Testing divand in 2 dimensions with independent verification.

using Base.Test

# grid of background field (its size should be odd)
xi,yi = ndgrid(linspace(0,1,15),linspace(0,1,15))

# mask (all points are valid)
mask = trues(xi)

# metric (inverse of the resolution)
pm = ones(xi) / (xi[2,1]-xi[1,1])
pn = ones(xi) / (yi[1,2]-yi[1,1])


# grid of observations
x = [0.5]
y = [0.5]
f = [1.]

# correlation length
len = 0.2

# normalized error variance
epsilon2 = 1.;


    fi,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);

	
        @test abs(mean(divand_diagHK(s))-0.5)<0.05
		@test abs(mean(divand_diagHKobs(s))-0.5)<0.05
		@test abs(mean(divand_diagHKobs(s,[1]))-0.5)<0.05
		@test abs(mean(divand_GCVKii(s))-0.5)<0.05
		@test abs(mean(divand_GCVKiiobs(s))-0.5)<0.05
		@test abs(mean(divand_residual(s,fi))-0.5)<0.05
		@test abs(mean(divand_residualobs(s,fi))-0.5)<0.05
		@test abs(mean(divand_erroratdatapoints(s))-0.5)<0.05
		@test abs(mean(divand_adaptedeps2(s,fi))-0.89)<0.05
		@test abs(divand_cvestimator(s,divand_residual(s,fi))-0.22)<0.05

# Copyright (C) 2014-2017 Alexander Barth <a.barth@ulg.ac.be>
#                         Jean-Marie Beckers <JM.Beckers@ulg.ac.be>
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
