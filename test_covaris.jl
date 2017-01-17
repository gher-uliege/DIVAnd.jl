using Base
using Base.Test
import Base

#type CovarIS{T} <: AbstractMatrix{eltype(T)}
#    IS::T
#end


Base.inv{T}(C::CovarIS{T}) = C.IS

Base.size{T}(C::CovarIS{T}) = size(IS)

Base.:*{T}(C::CovarIS{T},M) = C \ M


IS = [2 1; 1 2]

n = 2;

#A = sprandn(n,n,0.5);
#IS = A*A' + sparse_diag(n);
det(IS)

C = CovarIS(IS);
C2 = inv(IS);

iC = inv(C);
@test iC ≈ IS


b = randn(n,2);

a = C*b;
a2 = C2*b;

@test a ≈ a2


a = C\b;
a2 = C2\b;

@test a ≈ a2


C = factorize(C);

a = C*b;
a2 = C2*b;

@test a ≈ a2


a = C\b;
a2 = C2\b;

@test a,a2)


@test diag(C) ≈ diag(C2)


