# Testing DIVAnd in 2 dimensions with independent verification.

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

# grid of background field (its size should be odd)
mask,(pm,pn),(xi,yi) = DIVAnd_squaredom(
    2,range(0, stop = 1, length = 9))

# grid of observations
x = [0.5]
y = [0.5]
f = [1.]

# correlation length
len = 0.6

# normalized error variance
epsilon2 = 1.;

function DIVAnd_error(args...)
    f,s = DIVAndrun(args...)
    return statevector_unpack(s.sv,diag(s.P))[1]
end

function DIVAnd_almostexacterror(args...)
    err,bjmb,fa,sa = DIVAnd_aexerr(args...)
    return err
end

errormethods = [
                # consistent error (expensive)
                DIVAnd_error,
                # clever poor man's error
                DIVAnd_cpme,
                # almost exact error
                DIVAnd_almostexacterror
                ]


for errormethod in errormethods
    #@show errormethod
    err = errormethod(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2);

    errmin,minloc = findmin(err)

    if VERSION >= v"0.7.0-beta.0"
        for i = 1:ndims(mask)
            @test minloc[i] == (size(xi,i)+1) ÷ 2
        end
    else
        # should be the middle of the domain
        minsub = ind2sub(size(mask),minloc)

        for i = 1:ndims(mask)
            @test minsub[i] == (size(xi,i)+1) ÷ 2
        end
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
