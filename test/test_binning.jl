using Test

@test DIVAnd.unstagger([1, 2, 3]) == [0.5, 1.5, 2.5, 3.5]


@test DIVAnd.unstagger(1:3) == 0.5:1.0:3.5


x = 1:0.1:5

@test DIVAnd.findin(x, 2.1234) == 12
@test DIVAnd.findin(collect(x), 2.1234) == 12

