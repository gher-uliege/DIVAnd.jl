

trans,invtrans = Anam.loglin(10.)

x = linspace(0.03,200,100)
x2 = invtrans.(trans.(x))

@test x ≈ x2
