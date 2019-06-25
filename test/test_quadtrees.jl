if VERSION >= v"0.7.0-beta.0"
    using Test
    using Random
    using Dates
else
    using Base.Test
end

using DIVAnd
import DIVAnd.Quadtrees

# Quadtrees = DIVAnd.Quadtrees

if VERSION >= v"0.7.0-beta.0"
   Random.seed!(123)
else
   srand(123)
end

X = [0  0;
     1  0;
     1  1;
     0  1;
     0  0]

qt = DIVAnd.Quadtrees.QT(copy(X'),collect(1:size(X,1)))
attribs_res = DIVAnd.Quadtrees.within(qt,[0,0],[0.1,0.1])
@test attribs_res == [1,5]

@test DIVAnd.Quadtrees.within_count(qt,[0,0],[0.1,0.1]) == 2


DIVAnd.Quadtrees.within(qt,[0,0],[0.1,0.1])

@test [DIVAnd.Quadtrees.bitget(42,i) for i = 6:-1:1] == [true,false,true,false,true,false]

@test DIVAnd.Quadtrees.inside([0,0],[1,1],[0.5,0.5]) == true
@test DIVAnd.Quadtrees.inside([0,0],[1,1],[1.5,1.5]) == false


@test DIVAnd.Quadtrees.intersect([0,0],[1,1],[0.5,0.5],[2,2]) == true
@test DIVAnd.Quadtrees.intersect([0,0],[1,1],[1.5,1.5],[2,2]) == false
# one rectange contains the other
@test DIVAnd.Quadtrees.intersect([0,0],[1,1],[-1,-1],[2,2]) == true
# one very wide and very narrow rectangle
@test Quadtrees.intersect([-100,-1],[100,1],[-1,-100],[1,100]) == true

#@test_throws ArgumentError intersect([0,0],[1,1],[0.5,0.5],[2,2,3])

qt = DIVAnd.Quadtrees.QT(Int,[0.,0.],[1.,1.])
DIVAnd.Quadtrees.add!(qt,[0.1,0.1],1)
DIVAnd.Quadtrees.add!(qt,[0.2,0.2],2)
DIVAnd.Quadtrees.add!(qt,[0.7,0.7],3)
DIVAnd.Quadtrees.add!(qt,[0.9,0.1],4)

DIVAnd.Quadtrees.split!(qt)


@test DIVAnd.Quadtrees.isleaf(qt) == false
@test DIVAnd.Quadtrees.isleaf(qt.children[1]) == true

X = rand(2,10000)
attribs = collect(1:size(X,2))
qt2 = DIVAnd.Quadtrees.QT(X,attribs)
DIVAnd.Quadtrees.rsplit!(qt2,5)

#rplot(qt2)
#plot(X[1,:],X[2,:],"b.")

xmin = [0.3,0.3]
xmax = [0.51,0.51]

function simplesearch(X,attribs,xmin,xmax)
    sel = trues(size(X,2))
    for j = 1:size(X,1)
        sel = sel .& (xmin[j] .<= X[j,:] .<= xmax[j])
    end
    ind = findall(sel)
    return X[:,ind],attribs[ind]
end


xref,attribs_ref = simplesearch(X,attribs,xmin,xmax)

attribs_res = DIVAnd.Quadtrees.within(qt2,xmin,xmax)

#@show xref

@test sort(attribs_ref) == sort(attribs_res)


# progressively add all points
qt3 = DIVAnd.Quadtrees.QT(Int,[0.,0.],[1.,1.])

for i = 1:size(X,2)
    DIVAnd.Quadtrees.add!(qt3,X[:,i],i)
end

@test DIVAnd.Quadtrees.count(qt3) == size(X,2)

attribs_res = DIVAnd.Quadtrees.within(qt3,xmin,xmax)
@test sort(attribs_ref) == sort(attribs_res)


# Test in 1D - 4D

for n = 1:4
    #@show n
    Xtest = rand(n,100)
    attribs_ = collect(1:size(Xtest,2))

    qtND = DIVAnd.Quadtrees.QT(Xtest,attribs_)
    #@show qtND.points[1:5,:]

    DIVAnd.Quadtrees.rsplit!(qtND)

    @test DIVAnd.Quadtrees.ndims(qtND) == n
    @test DIVAnd.Quadtrees.count(qtND) == size(Xtest,2)

    xmin_ = fill(0.0,(n,))
    xmax_ = fill(0.5,(n,))

    xref_,attribs_ref_ = simplesearch(Xtest,attribs_,xmin_,xmax_)
    attribs_res_ = DIVAnd.Quadtrees.within(qtND,xmin_,xmax_)

    @test sort(attribs_ref_) == sort(attribs_res_)

    s = IOBuffer()
    show(s,qtND)
    @test occursin("Node",String(take!(s)))


end

# duplicates
x = Float64[i for i = 1:10, j = 1:11]
y = Float64[j for i = 1:10, j = 1:11]

value = x+y;
sel = [1:length(x);1];
dup = DIVAnd.Quadtrees.checkduplicates((x[sel],y[sel]),value[sel],[0.1,0.1],0.01)
@test dup == [[1,length(sel)]]

dup = DIVAnd.Quadtrees.checkduplicates((x[:],y[:]),value[:],
                                       (x[1:2],y[1:2]),value[1:2],[0.1,0.1],0.01)
@test dup == [[1],[2]]
