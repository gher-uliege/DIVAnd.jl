using Test
using LinearAlgebra

# Testing sparse and MatFun operators.

for opertype in [Val{:sparse}, Val{:MatFun}]



    field = randn(21, 30, 10)

    Sdiff = oper_diff(opertype, size(field), 1)
    f1 = field[2:end, :, :] - field[1:end-1, :, :]
    f2 = Sdiff * field[:]
    @test f1[:] ≈ f2

    Sdiff = oper_diff(opertype, size(field), 2)
    f1 = field[:, 2:end, :] - field[:, 1:end-1, :]
    f2 = Sdiff * field[:]
    @test f1[:] ≈ f2

    Sdiff = oper_diff(opertype, size(field), 3)
    f1 = field[:, :, 2:end] - field[:, :, 1:end-1]
    f2 = Sdiff * field[:]
    @test f1[:] ≈ f2

    # cyclic

    # dim 1 cyclic
    fc = vcat(field, reshape(field[1, :, :], (1, size(field, 2), size(field, 3))))
    dfc = fc[2:end, :, :] - fc[1:end-1, :, :]

    Sdiff = oper_diff(opertype, size(field), 1, true)
    f2 = Sdiff * field[:]
    @test dfc[:] ≈ f2


    # dim 2 cyclic
    fc = hcat(field, reshape(field[:, 1, :], (size(field, 1), 1, size(field, 3))))
    f1 = fc[:, 2:end, :] - fc[:, 1:end-1, :]

    Sdiff = oper_diff(opertype, size(field), 2, true)
    f2 = Sdiff * field[:]
    @test f1[:] ≈ f2

    # stagger

    Sdiff = oper_stagger(opertype, size(field), 1)
    f1 = (field[2:end, :, :] + field[1:end-1, :, :]) / 2
    f2 = Sdiff * field[:]
    @test f1[:] ≈ f2


    # dim 1 cyclic
    fc = vcat(field, reshape(field[1, :, :], (1, size(field, 2), size(field, 3))))
    f1 = (fc[2:end, :, :] + fc[1:end-1, :, :]) / 2

    Sdiff = oper_stagger(opertype, size(field), 1, true)
    f2 = Sdiff * field[:]
    @test f1[:] ≈ f2

    # shifting

    Sdiff = oper_shift(opertype, size(field), 1)
    f1 = field[2:end, :, :]
    f2 = Sdiff * field[:]
    @test f1[:] ≈ f2

    # trimming

    Sdiff = oper_trim(opertype, size(field), 1)
    f1 = field[2:end-1, :, :]
    f2 = Sdiff * field[:]
    @test f1[:] ≈ f2

    # sparse pack

    masktest = rand(Bool, size(field))
    masktest[1] = true
    f1 = field[masktest]
    f2 = oper_pack(opertype, masktest) * field[:]

    @test f1[:] ≈ f2


    # sparse interp

    masktest = trues(size(masktest))
    gridindices = [2.5; 2; 2]
    Hoperator, out, outbbox = sparse_interp(masktest, gridindices)
    f1 = (field[2, 2, 2] + field[3, 2, 2]) / 2
    f2 = Hoperator * field[:]
    @test [f1] ≈ f2

    # sparse interp

    gridindices = [2.5 3; 2 3; 2 3]
    Hoperator, out, outbbox = sparse_interp(masktest, gridindices)
    f1 = [(field[2, 2, 2] + field[3, 2, 2]) / 2; field[3, 3, 3]]
    f2 = Hoperator * field[:]

    @test f1 ≈ f2

    #  laplacian 1D

    coordx1 = 2 * collect(1:5)
    masktest = trues(size(coordx1))
    pm = ones(size(coordx1)) / 2

    DD = DIVAnd_laplacian(opertype, masktest, (pm,), ones(size(coordx1)), [false])
    field = 2 * coordx1.^2
    Df1 = 4
    Df2 = reshape(DD * field[:], size(masktest))
    Df2 = Df2[2:end-1]

    @test Df1 ≈ Df2[1]

    # value outside of grid
    Hoperator, out, outbbox = sparse_interp(trues(2, 2), [-10; 1])
    @test out == [true]

    # sparse gradient

    coordx1, coordx2 = ndgrid(2 * collect(1:4), 3 * collect(1:3))
    masktest = trues(size(coordx1))
    pm = ones(size(coordx1)) / 2
    pn = ones(size(coordx1)) / 3
    Dx, Dy = DIVAnd_gradient(opertype, masktest, (pm, pn))
    field = 2 * coordx1 + coordx2
    Df1 = 2 * ones(3, 3)
    Df2 = Dx * field[:]
    @test Df1[:] ≈ Df2

    # laplacian

    coordx1, coordx2 = ndgrid(2 * collect(1:4), 3 * collect(1:3))
    masktest = trues(size(coordx1))
    pm = ones(size(coordx1)) / 2
    pn = ones(size(coordx1)) / 3


    coordx1, coordx2 = ndgrid(collect(1:4), collect(1:3))
    masktest = trues(size(coordx1))
    pm = ones(size(coordx1))
    pn = ones(size(coordx1))
    DD = DIVAnd_laplacian(
        opertype,
        masktest,
        (pm, pn),
        ones(size(masktest)),
        [false, false],
    )
    field = 2 * coordx1.^2 + coordx2
    Df1 = 4.
    Df2 = reshape(DD * field[:], size(masktest))
    Df2 = Df2[2:end-1, 2:end-1]

    @test Df1 ≈ Df2[1]

end


# adjoint

operatortype = Val{:MatFun}
f = randn(10)
S = oper_diff(operatortype, size(f), 1)
a = randn(size(S, 1))
b = randn(size(S, 2))
@test a ⋅ (S * b) ≈ b ⋅ (S' * a)

sz = (10, 9, 8)
# diff and non cyclic
S = oper_diff(operatortype, sz, 3)
a = randn(size(S, 1))
b = randn(size(S, 2))
@test a ⋅ (S * b) ≈ b ⋅ (S' * a)

# diff and cyclic
S = oper_diff(operatortype, sz, 3, true)
a = randn(size(S, 1))
b = randn(size(S, 2))
@test a ⋅ (S * b) ≈ b ⋅ (S' * a)

# shift and non cyclic
S = oper_shift(operatortype, sz, 3)
a = randn(size(S, 1))
b = randn(size(S, 2))
@test a ⋅ (S * b) ≈ b ⋅ (S' * a)

# shift and cyclic
S = oper_shift(operatortype, sz, 3, true)
a = randn(size(S, 1))
b = randn(size(S, 2))
@test a ⋅ (S * b) ≈ b ⋅ (S' * a)

# stagger and non cyclic
S = oper_stagger(operatortype, sz, 3)
a = randn(size(S, 1))
b = randn(size(S, 2))
@test a ⋅ (S * b) ≈ b ⋅ (S' * a)

# stagger and cyclic
S = oper_stagger(operatortype, sz, 3, true)
a = randn(size(S, 1))
b = randn(size(S, 2))
@test a ⋅ (S * b) ≈ b ⋅ (S' * a)

# trim
S = oper_trim(operatortype, sz, 3)
a = randn(size(S, 1))
b = randn(size(S, 2))
@test a ⋅ (S * b) ≈ b ⋅ (S' * a)

# pack
mask = rand(10, 11) .> 0.5
S = oper_pack(operatortype, mask)
a = randn(size(S, 1))
b = randn(size(S, 2))
@test a ⋅ (S * b) ≈ b ⋅ (S' * a)


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
