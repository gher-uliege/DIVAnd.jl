# A simple example of DIVAnd in 1 dimensions
# with observations from an analytical function,

using DIVAnd
using Compat: @info, range
using PyPlot

include("./prep_dirs.jl")

# observations with points outside
x = collect(range(0,stop=1,length=7))
f = sin.(3*pi*x) ;

# final grid

xi=collect(range(-0.1,stop=1.1,length=220));

# reference field
fref = sin.(xi*6*pi) ;

# all points are valid points
mask = trues(size(xi));

pm = ones(size(xi)) / (xi[2]-xi[1]);

# correlation length
len = 0.05;

# obs. error variance normalized by the background error variance
epsilon2 = 1.;

m = Int(ceil(1+1/2))
# alpha is the (m+1)th row of the Pascal triangle:
# m=0         1
# m=1       1   1
# m=1     1   2   1
# m=2   1   3   3   1
# ...

alpha = [binomial(m,k) for k = 0:m];
@show(alpha)
# fi is the interpolated field

firef,s = DIVAndrun(mask,(pm,),(xi,),(x,),f,len,epsilon2;);

fi4,s = DIVAndrun(mask,(pm,),(xi,),(x,),f,len,epsilon2;alpha=alpha);

alpha = [binomial(m,k) for k = 0:m];
alpha = 2 * alpha
# fi is the interpolated field
fi1,s = DIVAndrun(mask,(pm,),(xi,),(x,),f,len,epsilon2;alpha=alpha);

alpha = [binomial(m,k) for k = 0:m];
alpha[1] = 0;
fi2,s = DIVAndrun(mask,(pm,),(xi,),(x,),f,len,epsilon2;alpha=alpha);

alpha = [binomial(m,k) for k = 0:m];
alpha[2] = 0;
fi3,s = DIVAndrun(mask,(pm,),(xi,),(x,),f,len,epsilon2;alpha=alpha);


plot(xi,fi1,".")
plot(x,f,"o",label="Observations")
plot(xi,fi2,"-")
plot(xi,fi3,":")
plot(xi,fi4,"+")
plot(xi,firef,"_",label="Reference")
legend()
figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => "_results.png")))
savefig(figname)
@info "Created figure " * figname


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
