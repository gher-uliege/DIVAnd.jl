# Testing DIVAnd in 2 dimensions with independent verification.

import DIVAnd
if VERSION >= v"0.7.0-beta.0"
    using Test
    using Statistics
else
    using Base.Test
end

# grid of background field (its size should be odd)
mask,(pm,pn),(xi,yi) = DIVAnd_squaredom(
    2,range(0.0,stop=1.0,length=15))

# grid of observations
x = [0.5]
y = [0.5]
f = [1.]

# correlation length
len = 0.2

# normalized error variance
epsilon2 = 1.;

fi,s = DIVAnd.DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;alphabc=0);

@test abs(mean(DIVAnd.DIVAnd_diagHK(s))-0.5)<0.05
@test abs(mean(DIVAnd.DIVAnd_diagHKobs(s))-0.5)<0.05
@test abs(mean(DIVAnd.DIVAnd_diagHKobs(s,[1]))-0.5)<0.05
@test abs(mean(DIVAnd.DIVAnd_GCVKii(s))-0.5)<0.05
@test abs(mean(DIVAnd.DIVAnd_GCVKiiobs(s))-0.5)<0.05
@test abs(mean(DIVAnd.DIVAnd_residual(s,fi))-0.5)<0.05
@test abs(mean(DIVAnd.DIVAnd_residualobs(s,fi))-0.5)<0.05
@test abs(mean(DIVAnd.DIVAnd_erroratdatapoints(s))-0.5)<0.05
@test abs(mean(DIVAnd.DIVAnd_adaptedeps2(s,fi))-0.89)<0.15
@test abs(DIVAnd.DIVAnd_cvestimator(s,DIVAnd.DIVAnd_residual(s,fi))-0.22)<0.05

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
