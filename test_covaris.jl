#module covarIS


using Base
using Base.Test
#export test


type CovarIS{T} <: AbstractMatrix{T}
    IS:: AbstractMatrix{T}
    factors
end

function CovarIS{T}(IS::AbstractMatrix{T})
    factors = nothing
    CovarIS(IS,factors)
end


Base.inv{T}(C::CovarIS{T}) = C.IS

Base.size{T}(C::CovarIS{T}) = size(IS)

function Base.:*{T}(C::CovarIS{T}, M::AbstractMatrix{Float64}) 
    if C.factors != nothing
        return C.factors \ M
    else
        return C.IS \ M
    end
end

function Base.:*{T}(C::CovarIS{T}, M::AbstractVector{Float64}) 
    if C.factors != nothing
        return C.factors \ M
    else
        return C.IS \ M
    end
end


function Base.getindex{T}(C::CovarIS{T}, i::Int,j::Int)
    ei = zeros(eltype(C.IS),size(C,1)); ei[i] = 1
    ej = zeros(eltype(C.IS),size(C,1)); ej[j] = 1

    return (ej'*(C*ei))[1]
end


Base.:\{T}(C::CovarIS{T}, M::AbstractArray{Float64,2}) = C.IS * M

function factorize!{T}(C::CovarIS{T})
    C.factors = cholfact(Symmetric(C.IS))
#    C.factors = cholfact(Symmetric(C.IS), Val{true})
#    C.factors = cholfact(Hermitian(C.IS), Val{true})
#    C.factors = cholfact(Hermitian(C.IS))
end

#function test()

IS = sparse([2. 0.1; 0.1 2.])
#IS = [2. 0; 0. 2.]

n = 2;

#A = sprandn(n,n,0.5);
#IS = A*A' + sparse_diag(n);
det(IS)

C = CovarIS(IS);
C2 = inv(full(IS));

iC = inv(C);
@test iC ≈ IS


b = randn(n,2);

a = C*b;
a2 = C2*b;

@test a ≈ a2

v = randn(n);

a = C*v;
a2 = C2*v;

@test a ≈ a2


a = C\b;
a2 = C2\b;

@test a ≈ a2


factorize!(C);

a = C*b;
a2 = C2*b;

@test a ≈ a2


a = C\b;
a2 = C2\b;

@test a ≈ a2

@test C[1,1] ≈ C2[1,1]

#@test diag(C) ≈ diag(C2)

@show typeof(C[1,1])
@which diag(C)

@show diag(C)


#end
#end
