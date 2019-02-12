# Testing DIVAnd in 2 dimensions with independent verification.

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
using DIVAnd

scalefactor = (3.,3)

# grid of background field

#mask,(pm,pn),(xi,yi) = DIVAnd_rectdom(collect(-10:8.),collect(0:3.))
mask,(pm,pn),(xi,yi) = DIVAnd_rectdom(collect(0:3.),collect(0:3.))
mask[:,end-1:end] .= false
mask[2,end-1] = true

#mask[:,1] .= false

epsilon = 1e-10;

# grid of observations

x,y = ndgrid([1., 3.],
             [0.99999])

v = copy(x)/maximum(x)
#v = ones(size(x))
#v[:,2] = reverse(v[:,2])

x = x[:]
y = y[:]
v = v[:]

lenx = 300.;
leny = 6.;

epsilon2 = 0.0001;

#,err,s
#va,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2,primal=true,alphabc=0,alpha=[1,3,3,1])
#va,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2,primal=true,alphabc=0,alpha=[0,1,1])
#va,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2,primal=true,alphabc=0,alpha=[0,1,0])
#va,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,(lenx,leny),epsilon2,primal=true,alphabc=0)


#va2,s2 = DIVAndrun(copy(mask'),(copy(pn'),copy(pm')),(copy(yi'),copy(xi')),(y,x),v,(leny,lenx),epsilon2,primal=true,alphabc=0)

len = DIVAnd.len_harmonize((lenx,leny),mask)
pmn = (pm,pn)





function derivative2!(dim,mask,pmn,len,va,D)
    pm,pn = pmn
    sz = size(mask)
    if dim == 1
        for j = 1:sz[2]
            for i = 2:sz[1]-1
                if mask[i-1,j] && mask[i,j] && mask[i+1,j]
                    D[i,j] += len[1][i,j]^2 * pm[i,j]^2 * (va[i-1,j] - 2*va[i,j] + va[i+1,j])
                end
            end
        end
    else
        for j = 2:sz[2]-1
            for i = 1:sz[1]
                if mask[i,j-1] && mask[i,j] && mask[i,j+1]
                    D[i,j] += len[2][i,j]^2 * pn[i,j]^2 * (va[i,j-1] - 2*va[i,j] + va[i,j+1])
                end
            end
        end
    end
    return D
end


derivative2(dim,mask,pmn,len,va) = derivative2!(dim,mask,pmn,len,va,zeros(size(mask)))

sz = (100,101)
mask = rand(sz...) .> 0.1
va = randn(sz...)
pmn = (2*ones(sz),3*ones(sz))
len = (0.2*ones(sz),0.4*ones(sz))

D1 = derivative2(1,mask,pmn,len,va)
D2 = derivative2(2,mask,pmn,len,va)

@test D1 ≈ DIVAnd.derivative2n(1,mask,pmn,len,va)
@test D2 ≈ DIVAnd.derivative2n(2,mask,pmn,len,va)

S1 = DIVAnd.sparse_derivative2n(1,mask,pmn,len)
@test D1 ≈ reshape(S1 * va[:],size(mask))

S2 = DIVAnd.sparse_derivative2n(2,mask,pmn,len)
@test D2 ≈ reshape(S2 * va[:],size(mask))

#nothing

# Copyright (C) 2014, 2016 Alexander Barth <a.barth@ulg.ac.be>
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
