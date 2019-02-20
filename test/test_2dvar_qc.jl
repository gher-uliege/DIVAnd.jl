# A simple example of DIVAnd in 2 dimensions
# with observations from an analytical function.
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
using DIVAnd


# observations
# same random set used as sometimes qc flags vary depending on actual noise
if VERSION >= v"0.7.0-beta.0"
   Random.seed!(11)
else
   srand(11)
end
x = rand(150);
y = rand(150);

# Put two points in specific locations

x[1]=0.25
y[1]=0.75

x[2]=0.75
y[2]=0.25


f = sin.(x*2*pi) .* sin.(y*2*pi);


f=f+0.25*randn(150);

# Now fake some mix up in  two points coordinates

x[2]=0.25
y[1]=0.75

x[1]=0.75
y[2]=0.25



# final grid
mask,(pm,pn),(xi,yi) = DIVAnd_squaredom(2,range(0,stop=1,length=20))

# correlation length
len = 0.1;

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

# fi is the interpolated field
fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),f,len,epsilon2;alphabc=0);


for method in [0, 1, 3, 4]
    qcval = DIVAnd_qc(fi,s,method)

    if method==4
        # Provide fake THETA value
        qcval=qcval*4
    end

    # Find suspect points
    sp=findall(x-> x.>9,qcval)
    @test sum(sp)==3
end

qcval_2 =
    @static if VERSION >= v"0.7.0"
        @test_logs (:warn,r".*not defined.*") match_mode=:any DIVAnd_qc(fi,s,2)
    else
        @test_warn r".*not defined.*" DIVAnd_qc(fi,s,2)
    end

@test all(qcval_2 .== 0)



# Copyright (C) 2014, 2017 Alexander Barth <a.barth@ulg.ac.be>
#                         Jean-Marie Beckers   <JM.Beckers@ulg.ac.be>
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
