if VERSION >= v"0.7.0-beta.0"
    using Test
    using LinearAlgebra
else
    using Base.Test
end

# Testing sparse and MatFun operators.

for operatortype in [Val{:sparse}, Val{:MatFun}]



    f = randn(21,30,10)

    Sdiff = oper_diff(operatortype,size(f),1)
    f1 = f[2:end,:,:] - f[1:end-1,:,:]
    f2 = Sdiff*f[:]
    @test f1[:] ≈ f2

    Sdiff = oper_diff(operatortype,size(f),2)
    f1 = f[:,2:end,:] - f[:,1:end-1,:]
    f2 = Sdiff*f[:]
    @test f1[:] ≈ f2

    Sdiff = oper_diff(operatortype,size(f),3)
    f1 = f[:,:,2:end] - f[:,:,1:end-1]
    f2 = Sdiff*f[:]
    @test f1[:] ≈ f2

    # cyclic

    # dim 1 cyclic
    fc =
        if VERSION >= v"0.7.0-beta.0"
            cat(f,reshape(f[1,:,:],(1,size(f,2),size(f,3))), dims=1)
        else            
            cat(1,f,reshape(f[1,:,:],(1,size(f,2),size(f,3))))
        end
    dfc = fc[2:end,:,:] - fc[1:end-1,:,:]

    Sdiff = oper_diff(operatortype,size(f),1,true)
    f2 = Sdiff*f[:]
    @test dfc[:] ≈ f2


    # dim 2 cyclic
    fc =
        if VERSION >= v"0.7.0-beta.0"
            cat(f,reshape(f[:,1,:],(size(f,1),1,size(f,3))), dims = 2)
        else
            cat(2,f,reshape(f[:,1,:],(size(f,1),1,size(f,3))))
        end
    f1 = fc[:,2:end,:] - fc[:,1:end-1,:]

    Sdiff = oper_diff(operatortype,size(f),2,true)
    f2 = Sdiff*f[:]
    @test f1[:] ≈ f2

    # stagger

    Sdiff = oper_stagger(operatortype,size(f),1)
    f1 = (f[2:end,:,:] + f[1:end-1,:,:])/2
    f2 = Sdiff*f[:]
    @test f1[:] ≈ f2


    # dim 1 cyclic
    fc =
        if VERSION >= v"0.7.0-beta.0"            
            cat(f,reshape(f[1,:,:],(1,size(f,2),size(f,3))), dims=1)
        else
            cat(1,f,reshape(f[1,:,:],(1,size(f,2),size(f,3))))
        end
    f1 = (fc[2:end,:,:] + fc[1:end-1,:,:])/2

    Sdiff = oper_stagger(operatortype,size(f),1,true)
    f2 = Sdiff*f[:]
    @test f1[:] ≈ f2

    # shifting

    Sdiff = oper_shift(operatortype,size(f),1)
    f1 = f[2:end,:,:]
    f2 = Sdiff*f[:]
    @test f1[:] ≈ f2

    # trimming

    Sdiff = oper_trim(operatortype,size(f),1)
    f1 = f[2:end-1,:,:]
    f2 = Sdiff*f[:]
    @test f1[:] ≈ f2

    # sparse pack

    mask = rand(Bool,size(f))
    mask[1] = true
    f1 = f[mask]
    f2 = oper_pack(operatortype,mask) * f[:]

    @test f1[:] ≈ f2


    # sparse interp

    mask = trues(size(mask))
    gridindices = [2.5; 2; 2]
    H,out,outbbox = sparse_interp(mask,gridindices)
    f1 = (f[2,2,2] + f[3,2,2])/2
    f2 = H*f[:]
    @test [f1] ≈ f2

    # sparse interp

    gridindices = [2.5 3; 2 3; 2 3]
    H,out,outbbox = sparse_interp(mask,gridindices)
    f1 = [ (f[2,2,2] + f[3,2,2])/2 ;  f[3,3,3] ]
    f2 = H*f[:]

    @test f1 ≈ f2

    #  laplacian 1D

    x1 = 2 * collect(1:5);
    mask = trues(size(x1));
    pm = ones(size(x1))/2;

    DD = DIVAnd_laplacian(operatortype,mask,(pm,),ones(size(x1)),[false]);
    f = 2*x1.^2;
    Df1 = 4;
    Df2 = reshape(DD * f[:], size(mask));
    Df2 = Df2[2:end-1]

    @test Df1 ≈ Df2[1]

    # value outside of grid
    H,out,outbbox = sparse_interp(trues(2,2),[-10; 1])
    @test out == [true]

    # sparse gradient

    x1,x2 = ndgrid(2*collect(1:4),3*collect(1:3))
    mask = trues(size(x1))
    pm = ones(size(x1))/2
    pn = ones(size(x1))/3
    Dx,Dy = sparse_gradient(operatortype,mask,(pm,pn))
    f = 2*x1 + x2
    Df1 = 2 * ones(3,3)
    Df2 = Dx * f[:]
    @test Df1[:] ≈ Df2

    # laplacian

    x1,x2 = ndgrid(2*collect(1:4),3*collect(1:3))
    mask = trues(size(x1))
    pm = ones(size(x1))/2
    pn = ones(size(x1))/3


    x1,x2 = ndgrid(collect(1:4),collect(1:3))
    mask = trues(size(x1))
    pm = ones(size(x1))
    pn = ones(size(x1))
    DD = DIVAnd_laplacian(operatortype,mask,(pm,pn),ones(size(mask)),[false,false])
    f = 2*x1.^2 + x2
    Df1 = 4.
    Df2 = reshape(DD * f[:], size(mask))
    Df2 = Df2[2:end-1,2:end-1]

    @test Df1 ≈ Df2[1]

end


# adjoint

operatortype = Val{:MatFun}
f = randn(10)
S = oper_diff(operatortype,size(f),1)
a = randn(size(S,1))
b = randn(size(S,2))
@test a ⋅ (S*b) ≈ b ⋅ (S'*a)

sz = (10,9,8)
# diff and non cyclic
S = oper_diff(operatortype,sz,3)
a = randn(size(S,1))
b = randn(size(S,2))
@test a ⋅ (S*b) ≈ b ⋅ (S'*a)

# diff and cyclic
S = oper_diff(operatortype,sz,3,true)
a = randn(size(S,1))
b = randn(size(S,2))
@test a ⋅ (S*b) ≈ b ⋅ (S'*a)

# shift and non cyclic
S = oper_shift(operatortype,sz,3)
a = randn(size(S,1))
b = randn(size(S,2))
@test a ⋅ (S*b) ≈ b ⋅ (S'*a)

# shift and cyclic
S = oper_shift(operatortype,sz,3,true)
a = randn(size(S,1))
b = randn(size(S,2))
@test a ⋅ (S*b) ≈ b ⋅ (S'*a)

# stagger and non cyclic
S = oper_stagger(operatortype,sz,3)
a = randn(size(S,1))
b = randn(size(S,2))
@test a ⋅ (S*b) ≈ b ⋅ (S'*a)

# stagger and cyclic
S = oper_stagger(operatortype,sz,3,true)
a = randn(size(S,1))
b = randn(size(S,2))
@test a ⋅ (S*b) ≈ b ⋅ (S'*a)

# trim
S = oper_trim(operatortype,sz,3)
a = randn(size(S,1))
b = randn(size(S,2))
@test a ⋅ (S*b) ≈ b ⋅ (S'*a)

# pack
mask = rand(10,11) .> 0.5
S = oper_pack(operatortype,mask)
a = randn(size(S,1))
b = randn(size(S,2))
@test a ⋅ (S*b) ≈ b ⋅ (S'*a)


# Copyright (C) 2014,2016 Alexander Barth <a.barth@ulg.ac.be>
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
