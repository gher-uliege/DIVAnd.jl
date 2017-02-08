# Testing divand in 2 dimensions with independent verification.

using Base.Test

# grid of background field (its size should be odd)
xi,yi = ndgrid(linspace(0,1,9),linspace(0,1,9))

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
len = .15;

# normalized error variance
epsilon2 = 1.;

function divand_error(args...)
    f,s = divandrun(args...)
    return statevector_unpack(s.sv,diag(s.P))[1]
end

function divand_almostexacterror(args...)
    err,bjmb,fa,sa = divand_aexerr(args...)
    return err
end

errormethods = [
                # consistent error (expensive)
                divand_error,
                # clever poor man's error
                divand_cpme,
                # almost exact error
                #divand_almostexacterror
                ]


for errormethod in errormethods
    #@show errormethod
    err = errormethod(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);

    errmin,minloc = findmin(err)

    # should be the middle of the domain
    minsub = ind2sub(size(mask),minloc)

    for i = 1:ndims(mask)
        @test minsub[i] == (size(xi,i)+1) ÷ 2
    end
end

# Copyright (C) 2014-2017 Alexander Barth <a.barth@ulg.ac.be>
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
