using Base.Test

# Testing sparse operators.

f = randn(21,30,10)

S = sparse_diff(size(f),1)
f1 = f[2:end,:,:] - f[1:end-1,:,:]
f2 = S*f[:]
@test f1[:] ≈ f2

S = sparse_diff(size(f),2)
f1 = f[:,2:end,:] - f[:,1:end-1,:]
f2 = S*f[:]
@test f1[:] ≈ f2

S = sparse_diff(size(f),3)
f1 = f[:,:,2:end] - f[:,:,1:end-1]
f2 = S*f[:]
@test f1[:] ≈ f2

# cyclic

# dim 1 cyclic
fc = cat(1,f,reshape(f[1,:,:],(1,size(f,2),size(f,3))))
dfc = fc[2:end,:,:] - fc[1:end-1,:,:]

S = sparse_diff(size(f),1,true)
f2 = S*f[:]
@test dfc[:] ≈ f2


# dim 2 cyclic
fc = cat(2,f,reshape(f[:,1,:],(size(f,1),1,size(f,3))))
f1 = fc[:,2:end,:] - fc[:,1:end-1,:]

S = sparse_diff(size(f),2,true)
f2 = S*f[:]
@test f1[:] ≈ f2

# stagger

S = sparse_stagger(size(f),1)
f1 = (f[2:end,:,:] + f[1:end-1,:,:])/2
f2 = S*f[:]
@test f1[:] ≈ f2


# # dim 1 cyclic
# fc = cat(1,f,reshape(f[1,:,:],(1,size(f,2),size(f,3))))
# f1 = (fc[2:end,:,:] + fc[1:end-1,:,:])/2

# S = sparse_stagger(size(f),1,true)
# f2 = S*f[:]
# @test f1[:] ≈ f2

# # shifting

# S = sparse_shift(size(f),1)
# f1 = f[2:end,:,:]
# f2 = S*f[:]
# @test f1[:] ≈ f2

# # trimming

# S = sparse_trim(size(f),1)
# f1 = f[2:end-1,:,:]
# f2 = S*f[:]
# @test f1[:] ≈ f2

# # sparse pack

# mask = rand(size(f)) .> 0
# f1 = f[mask]
# f2 = sparse_pack(mask) * f[:]

# @test f1[:] ≈ f2


# # sparse interp

# mask = trues(size(mask))
# I = [2.5 2 2]'
# H,out,outbbox = sparse_interp(mask,I)
# f1 = (f[2,2,2] + f[3,2,2])/2
# f2 = H*f[:]
# @test [f1] ≈ f2

# # sparse interp

# I = [2.5 3; 2 3; 2 3]
# H,out,outbbox = sparse_interp(mask,I)
# f1 = [ (f[2,2,2] + f[3,2,2])/2 ;  f[3,3,3] ]
# f2 = H*f[:]

# @test f1 ≈ f2

# #  laplacian 1D

# x1 = 2 * collect(1:5);
# mask = trues(size(x1));
# pm = ones(size(x1))/2;

# DD = divand_laplacian(mask,(pm,),ones(size(x1)),[false]);
# f = 2*x1.^2;
# Df1 = 4;
# Df2 = reshape(DD * f[:], size(mask));
# Df2 = Df2[2:end-1]

# @test Df1 ≈ Df2[1]


# # sparse gradient

# x1,x2 = ndgrid(2*collect(1:4),3*collect(1:3))
# mask = trues(size(x1))
# pm = ones(size(x1))/2
# pn = ones(size(x1))/3
# Dx,Dy = sparse_gradient(mask,(pm,pn))
# f = 2*x1 + x2
# Df1 = 2 * ones(3,3)
# Df2 = Dx * f[:]
# @test Df1[:] ≈ Df2

# # laplacian

# x1,x2 = ndgrid(2*collect(1:4),3*collect(1:3))
# mask = trues(size(x1))
# pm = ones(size(x1))/2
# pn = ones(size(x1))/3


# x1,x2 = ndgrid(collect(1:4),collect(1:3))
# mask = trues(size(x1))
# pm = ones(size(x1))
# pn = ones(size(x1))
# DD = divand_laplacian(mask,(pm,pn),ones(size(mask)),[false,false])
# f = 2*x1.^2 + x2
# Df1 = 4.
# Df2 = reshape(DD * f[:], size(mask))
# Df2 = Df2[2:end-1,2:end-1]

# @test Df1 ≈ Df2[1]


# # Copyright (C) 2014,2016 Alexander Barth <a.barth@ulg.ac.be>
# #
# # This program is free software; you can redistribute it and/or modify it under
# # the terms of the GNU General Public License as published by the Free Software
# # Foundation; either version 2 of the License, or (at your option) any later
# # version.
# #
# # This program is distributed in the hope that it will be useful, but WITHOUT
# # ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# # FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# # details.
# #
# # You should have received a copy of the GNU General Public License along with
# # this program; if not, see <http://www.gnu.org/licenses/>.
