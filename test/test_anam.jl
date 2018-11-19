if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
include("../src/anamorphosis.jl")

# log/exp transform (with linear extension)

threshold = 10.
trans,invtrans = Anam.loglin(threshold)

x = range(0.03,stop = 200,length = 100)
x2 = invtrans.(trans.(x))

@test x ≈ x2
@test trans(threshold-10*eps(1.)) ≈ trans(threshold+10*eps(1.)) atol=100*eps(1.)

# with offset

trans,invtrans = Anam.loglin(threshold; epsilon = 0.05)

x = range(-0.03,stop = 200,length = 100)
x2 = invtrans.(trans.(x))

@test x ≈ x2
@test trans(threshold-10*eps(1.)) ≈ trans(threshold+10*eps(1.)) atol=100*eps(1.)

# logit transform

trans,invtrans = Anam.logit()

x = range(0.01, stop = 0.99, length = 100)
x2 = invtrans.(trans.(x))

@test x ≈ x2

# with custom range
trans,invtrans = Anam.logit(min = -2.1, max = 2.1)

@test 1 ≈ invtrans(trans(1))

x = range(-2, stop = 2,length = 100)
x2 = invtrans.(trans.(x))
@test x ≈ x2

# no transformation (default)

trans,invtrans = Anam.notransform()

x = range(0.03,stop = 200,length = 100)

@test collect(x) == invtrans.(trans.(x))
@test collect(x) == trans.(x)
