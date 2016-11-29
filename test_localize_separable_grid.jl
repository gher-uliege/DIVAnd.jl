using Base.Test


I = localize_separable_grid(([4],),ones(10),(2*collect(1:10),))
@show I
@test I[1] ≈ 2.


x1,y1 = ndgrid(2 * collect(1:5),collect(1:6))
x = (x1,y1)
xi = (3,3)
mask = trues(size(x1));

I = localize_separable_grid(xi,mask,x)

@test I ≈ [1.5; 3]
