# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

# true error variance of observation
epsilon2_true = 1.

if VERSION >= v"0.7.0-beta.0"
   Random.seed!(1234)
else
   srand(1234)
end
# observations
nobs = 99
x = rand(nobs);
y = rand(nobs);

# true length-scale
len_true = 0.5
f = sin.(π * x/len_true) .* cos.(π * y/len_true);

# add noise
f = f+sqrt(epsilon2_true) * randn(nobs);

# final grid
mask,(pm,pn),(xi,yi) = DIVAnd_rectdom(
    range(0,stop=1,length=14),range(0,stop=1,length=13))

# correlation length (first guess)
len = 0.1;

# obs. error variance normalized by the background error variance (first guess)
epsilon2 = 2.

# loop over all methods
for imeth = 0:3
    bestfactorl,bestfactore, cvval,cvvalues, x2Ddata,y2Ddata,cvinter,xi2D,yi2D = DIVAnd_cv(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2,2,3,imeth;alphabc=0);

    if VERSION >= v"0.7.0-beta.0"
        @test 0.5 < bestfactore*epsilon2/epsilon2_true <= 2
    else
        @test 0.5 < bestfactore*epsilon2/epsilon2_true < 2
    end
    @test 0.3 < bestfactorl*len/len_true < 3

    #@show bestfactore*epsilon2
    #@show bestfactorl*len
end

for imeth = 0:3
    bestfactor, cvval,cvvalues, logfactorse,cvinter,epsilon2inter = DIVAnd_cv(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2,0,3,imeth;alphabc=0);
    @test 1. < bestfactor < 1.3
    #@test 0.3 < bestfactorl*len/len_true < 3

    #@show bestfactor

end


for imeth = 0:3
    bestfactor, cvval,cvvalues, logfactorse,cvinter,epsilon2inter = DIVAnd_cv(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2,2,0,imeth;alphabc=0);
    #@test 0.5 < bestfactore*epsilon2/epsilon2_true < 2
    #@test 0.3 < bestfactorl*len/len_true < 3

    if VERSION >= v"0.7.0-beta.0"
        @test 1.5 < bestfactor < 1.8
    else
        @test 1.6 < bestfactor < 1.8
    end

    #@show bestfactor
    #@show bestfactore*epsilon2
    #@show bestfactorl*len
end

for imeth = 0:3
    bestfactor,cvvalues =
        @static if VERSION >= v"0.7.0"
            @test_logs (:warn,r".*no parameter optimisation.*") match_mode=:any DIVAnd_cv(
                mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2,0,0,imeth;alphabc=0);
        else
            @test_warn r".*no parameter optimisation.*" DIVAnd_cv(
                mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2,0,0,imeth;alphabc=0);
        end

    #@show bestfactor
    if VERSION >= v"0.7.0-beta.0"
        @test 0.8 < bestfactor < 1.1
    else
        @test 0.8 < bestfactor < 1.
    end
    #@show bestfactore*epsilon2
    #@show bestfactorl*len
end


# Copyright (C) 2014, 2017, 2018
#   Alexander Barth <a.barth@ulg.ac.be>
#   Jean-Marie Beckers <JM.Beckers@ulg.ac.be>
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
