using Base.Cartesian
if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
import DIVAnd

for N = 1:6
    @show N
    @eval begin

        function DIVAnd_laplacian_prepare2(mask::BitArray{$N},
                                           pmn::NTuple{$N,Vector{T}},
                                           nu::NTuple{$N,Vector{T}}) where T
            sz = size(mask)
            ivol = ones(T,sz)

            @nloops $N i d->1:sz[d] begin
                (@nref $N ivol i) = 1.
                @nexprs $N m->begin
                    (@nref $N ivol i) *= pmn[m][i_m]
                end
            end

            nus = ntuple(i -> zeros(T,sz),$N)::NTuple{$N,Array{T,$N}}

            # This heavily uses macros to generate fast code
            # In e.g. 3 dimensions
            # (@nref $N tmp i) corresponds to tmp[i_1,i_2,i_3]
            # (@nref $N nu_i l->(l==j?i_l+1:i_l)  corresponds to nu_i[i_1+1,i_2,i_3] if j==1

            # loop over all dimensions to create
            # nus[1] (nu stagger in the 1st dimension)
            # nus[2] (nu stagger in the 2nd dimension)
            # ...

            @nexprs $N j->begin
                tmp = nus[j]

                # loop over all spatio-temporal dimensions
                @nloops $N i k->(k == j ? (1:sz[k]-1) : (1:sz[k]))  begin
                    nu_i = nu[j]
                    # stagger nu
                    # e.g. 0.5 * (nu[1][2:end,:] + nu[1][1:end-1,:])
                    (@nref $N tmp i) = 0.5 * (nu_i[i_j+1] + nu_i[i_j])
                    # e.g. (pmn[2][2:end,:]+pmn[2][1:end-1,:]) ./ (pmn[1][2:end,:]+pmn[1][1:end-1,:])
                    @nexprs $N m->begin
                        pm_i = pmn[m]
                        if (m .== j)
                            (@nref $N tmp i) *= 0.5 * (pm_i[i_j+1] + pm_i[i_j])
                        else
                            (@nref $N tmp i) /= pm_i[i_m]
                        end

                        if !(@nref $N mask i) || !(@nref $N mask l->(l==j ? i_l+1 : i_l))
                            (@nref $N tmp i) = 0
                        end
                    end
                end
            end


            # for i = 1:length(sz)
            #     for j = 1:sz[i]-1
            #         nus[i][j] = 0.25 * (nu[i][j+1] + nu[i][j]) * (pmn[i][j+1] + pmn[i][j])


            #         for l = 1:length(sz)
            #             if l != i
            #                 nus[i][j] /= 0.5 * (pmn[l][j+1] + pmn[l][j])
            #             end
            #         end
            #     end

            # end

            return ivol,nus
        end


        function DIVAnd_laplacian_apply2!(mask,ivol,nus,x::AbstractArray{T,$N},Lx::AbstractArray{T,$N}) where T
            sz = size(x)
            Lx[:] = 0

            @inbounds @nloops $N i d->1:sz[d] begin
                (@nref $N Lx i) = T(0)
                @nexprs $N d1->begin
                    tmp2 = nus[d1]

                    if (@nref $N mask i)
                        if i_d1 < sz[d1]
                            @inbounds if (@nref $N mask d2->(d2==d1 ? i_d2+1 : i_d2))
                                @inbounds (@nref $N Lx i) += tmp2[i_d1] * ((@nref $N x d2->(d2==d1 ? i_d2+1 : i_d2)) - (@nref $N x i))
                            end
                        end

                        if i_d1 > 1
                            @inbounds if (@nref 2 mask d2->(d2==d1 ? i_d2-1 : i_d2))
                                @inbounds (@nref $N Lx i) -= tmp2[i_d1-1] * ((@nref $N x i) - (@nref $N x d2->(d2==d1 ? i_d2-1 : i_d2)))
                            end
                        end
                    end
                end

                (@nref $N Lx i) *= (@nref $N ivol i)
            end
        end



        function DIVAnd_laplacian_apply3!(mask,pmn,nu,x::AbstractArray{T,$N},Lx::AbstractArray{T,$N}) where T
            sz = size(x)

            # This heavily uses macros to generate fast code
            # In e.g. 3 dimensions
            # (@nref $N tmp i) corresponds to tmp[i_1,i_2,i_3]
            # (@nref $N nu_i l->(l==j?i_l+1:i_l)  corresponds to nu_i[i_1+1,i_2,i_3] if j==1

            # loop over all dimensions
            @inbounds @nloops $N i d->1:sz[d]  begin

                # initialize to zero Lx
                (@nref $N Lx i) = 0

                if (@nref $N mask i)
                    # for every dimension add the fluxes along the d1-th dimension
                    @nexprs $N d1->begin
                        # temporary array to reduce the number of indices in the following
                        nu_i = nu[d1]

                        #
                        if i_d1 < sz[d1]
                            if (@nref $N mask l->(l==d1 ? i_l+1 : i_l))

                                mytmp = 0.5 * (nu_i[i_d1+1] + nu_i[i_d1])
                                @nexprs $N m->begin
                                    pm_i = pmn[m]
                                    if (m .== d1)
                                        mytmp *= 0.5 * (pm_i[i_d1+1] + pm_i[i_d1])
                                    else
                                        mytmp /= pm_i[i_m]
                                    end
                                end

                                @inbounds (@nref $N Lx i) += mytmp * ((@nref $N x d2->(d2==d1 ? i_d2+1 : i_d2)) - (@nref $N x i))
                            end
                        end

                        if i_d1 > 1
                            if (@nref $N mask l->(l==d1 ? i_l-1 : i_l))

                                mytmp = 0.5 * (nu_i[i_d1] + nu_i[i_d1-1])
                                @nexprs $N m->begin
                                    pm_i = pmn[m]
                                    if (m .== d1)
                                        mytmp *= 0.5 * (pm_i[i_d1] + pm_i[i_d1-1])
                                    else
                                        mytmp /= pm_i[i_m]
                                    end
                                end


                                @inbounds (@nref $N Lx i) -= mytmp * ((@nref $N x i) - (@nref $N x d2->(d2==d1 ? i_d2-1 : i_d2)))
                            end
                        end
                    end

                    @nexprs $N m->begin
                        (@nref $N Lx i) *= pmn[m][i_m]
                    end
                end
            end
        end

end # @eval
end # for
function test_sp(mask,pmn,nu,x0,Nmax=1)
    sz = size(x0)
    x = x0[mask]

    L = DIVAnd.DIVAnd_laplacian(Val{:sparse},mask,pmn,nu,falses(ndims(mask)))

    @time @inbounds for nt = 1:Nmax
        x = L*x
    end

    Lx = zeros(x0)
    Lx[mask] = x
    return Lx
end


function test_sp_inplace(mask,pmn,nu,x0,Nmax=1)
    sz = size(x0)
    x = x0[mask]

    L = DIVAnd.DIVAnd_laplacian(Val{:sparse},mask,pmn,nu,falses(ndims(mask)))
    Lx = similar(x)

    @time @inbounds for nt = 1:Nmax
        #x = L*x
        A_mul_B!(Lx,L,x)
        x[:] = Lx
    end

    Lx2 = zeros(x0)
    Lx2[mask] = x
    return Lx2
end

function test_lap8(mask,pmn,nu,x0,Nmax=1)
    ivol,nus = DIVAnd.DIVAnd_laplacian_prepare(mask,pmn,nu)
    x = copy(x0)
    Lx = similar(x)

    @time @inbounds for nt = 1:Nmax
        DIVAnd.DIVAnd_laplacian_apply!(ivol,nus,x,Lx)
        x[:] = Lx
    end

    return x
end

function test_lapv(mask,pmn,nu,x0,Nmax=1)
    ivol,nus = DIVAnd_laplacian_prepare2(mask,pmn,nu)
    x = copy(x0)
    Lx = similar(x)

    @time @inbounds for nt = 1:Nmax
        #DIVAnd_laplacian_apply2!(mask,ivol,nus,x,Lx)
        DIVAnd.DIVAnd_laplacian_apply!(ivol,nus,x,Lx)
        x[:] = Lx
    end

    return x
end



function test_lapv3(mask,pmn,nu,x0,Nmax=1)
    x = copy(x0)
    Lx = similar(x)

    @time @inbounds for nt = 1:Nmax
        DIVAnd_laplacian_apply3!(mask,pmn,nu,x,Lx)
        x[:] = Lx
    end

    return x
end

#sz = (1000,1000)  # too large for a 32-bit OS
Nmax = 2

#for sz in [(20,),(100,100),(20,20,20),(5,5,5,5)]
sz = (200,200,200,5)
#sz = (20,20,20,4)
#sz = (5,5)
sz = (50,50,50)

x = randn(sz)
mask = trues(sz)

ij = DIVAnd.ndgrid([Float64.(1:s) for s in sz]...)
x = ij[1].^2

if length(sz) == 20
    pmn = (ij[1],ij[2]+ij[1]/10)
    nu = (ij[1]+2 * ij[2],ij[2]+3 * ij[1])

    mask[3:4,3:4] = false
else
    pmnv = ntuple(i -> (i+2.) * collect(1:sz[i]),length(sz))
    nuv = ntuple(i -> (i+2.) * collect(1:sz[i]).^2,length(sz))

    pmn = DIVAnd.ndgrid(pmnv...)
    nu = DIVAnd.ndgrid(nuv...)

end



Lx2v = test_lapv3(mask,pmnv,nuv,x,Nmax);

Lxsp0 = test_sp(mask,pmn,nu,x,Nmax);
include("../../src/override_ssmult.jl")
Lxsp = test_sp(mask,pmn,nu,x,Nmax);
Lx2 = test_lap8(mask,pmn,nu,x,Nmax);
Lxsp2 = test_sp_inplace(mask,pmn,nu,x,Nmax);


@test Lxsp ≈ Lx2 atol=1e-4
@test Lxsp ≈ Lx2v atol=1e-4

@test Lx2 ≈ Lx2v atol=1e-4

nothing
