using Test
using DIVAnd

gridindices = localize_separable_grid(([4],), trues((10,)), (2 * collect(1:10),))
@test gridindices[1] ≈ 2.

# 2D with one point

x1, y1 = ndgrid(2 * collect(1:5), collect(1:6))
x = (x1, y1)
xi = ([3], [3])
mask = trues(size(x1))

gridindices = localize_separable_grid(xi, mask, x)

@test gridindices ≈ [1.5; 3]


# 2D with 2 points

x1, y1 = ndgrid(2 * collect(1:5), collect(1:6))
x = (x1, y1)
xi = ([3, 4], [3, 5])
mask = trues(size(x1))

gridindices = localize_separable_grid(xi, mask, x)

@test gridindices ≈ [1.5 2.0; 3. 5.]


# 2D with 1 point outside

x1, y1 = ndgrid(range(0.5, stop = 1, length = 50),
                range(0., stop = 1, length = 30));
x = (x1, y1)
xi = ([0.2], [0.5])
mask = trues(size(x1))
gridindices = localize_separable_grid(xi, mask, x)
@test gridindices[1] < 1


x1, y1 = ndgrid(0.:10,0.:10)
x = (x1, y1)
iscyclic = (true,false)
xi = ([10.5, 5.], [2.,20000])
mask = trues(size(x1))
gridindices = localize_separable_grid(xi, mask, x, iscyclic)
@test gridindices[1,1] ≈ 11.5
@test gridindices[2,1] ≈ 3

@test gridindices[1,2] == -1
@test gridindices[2,2] == -1


#=
# benchmark

x1,y1,z1 = ndgrid(
    range(0,stop = 1, length = 100),
    range(0,stop = 1, length = 100),
    range(0,stop = 1, length = 10))

x = (x1,y1,z1)
m = 1_000_000

xi = (rand(Float64,m),rand(Float64,m),rand(Float64,m))
mask = trues(size(x1));

#@code_warntype localize_separable_grid(xi,mask,x)
using BenchmarkTools
gridindices = @btime localize_separable_grid(xi,mask,x);
# 3.840 s (38999071 allocations: 987.96 MiB)
# 692.904 ms (151 allocations: 72.44 MiB) -> after optimization
=#
