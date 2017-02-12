using Base.Test

function testprod(C,C2)

    n = size(C,2)

    # C times a matrix
    b = randn(n,2);
    a = C*b;
    a2 = C2*b;
    @test a ≈ a2

    # C times a matrix tranposed
    b = randn(2,n);
    a = C*b.';
    a2 = C2*b.';
    @test a ≈ a2

    # C times a matrix conjugate tranposed
    b = randn(2,n);
    a = C*b';
    a2 = C2*b';
    @test a ≈ a2

    # C times a vector
    v = randn(n);
    a = C*v;
    a2 = C2*v;
    @test a ≈ a2


    # diagonal
    @test diag(C) ≈ diag(C2)
end


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

testprod(C,C2)

# inverse of C times a matrix
b = randn(n,2);
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




# MatFun

M = randn(10,10)
MF = MatFun(size(M),x -> M*x,x -> M'*x)

x = randn(size(M,2))
A = randn(size(M,2),3)

@test size(M) == size(MF)
@test M*x ≈ MF*x
@test M'*x ≈ MF'*x
@test M*A ≈ MF*A



# CovarHPHt


P = randn(10,10); P = P*P';
H = randn(10,10)

A = CovarHPHt(P,H)
A2 = H*P*H'

testprod(A,A2)
