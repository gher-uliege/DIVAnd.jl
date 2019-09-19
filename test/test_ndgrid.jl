# test if ndgrid accepts different types of arguments

x, y = DIVAnd.ndgrid([1, 2, 3], [2, 3])
@test x[2, 1] == 2
@test y[2, 1] == 2


x2, y2 = DIVAnd.ndgrid([1, 2, 3], [2., 3.])
@test x == x2
@test y == y2

x3, y3 = DIVAnd.ndgrid(1:3, [2., 3.])
@test x == x2
@test y == y2
