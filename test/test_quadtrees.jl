if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end

import divand.Quadtrees

# Quadtrees = divand.Quadtrees

srand(123)

X = [0  0;
     1  0;
     1  1;
     0  1;
     0  0]

qt = divand.Quadtrees.QT(copy(X'),collect(1:size(X,1)))
attribs_res = divand.Quadtrees.within(qt,[0,0],[0.1,0.1])
@test attribs_res == [1,5]


@test [divand.Quadtrees.bitget(42,i) for i = 6:-1:1] == [true,false,true,false,true,false]

@test divand.Quadtrees.inside([0,0],[1,1],[0.5,0.5]) == true
@test divand.Quadtrees.inside([0,0],[1,1],[1.5,1.5]) == false


@test divand.Quadtrees.intersect([0,0],[1,1],[0.5,0.5],[2,2]) == true
@test divand.Quadtrees.intersect([0,0],[1,1],[1.5,1.5],[2,2]) == false
# one rectange contains the other
@test divand.Quadtrees.intersect([0,0],[1,1],[-1,-1],[2,2]) == true
#@test_throws ArgumentError intersect([0,0],[1,1],[0.5,0.5],[2,2,3])

qt = divand.Quadtrees.QT(Int,[0.,0.],[1.,1.])
divand.Quadtrees.add!(qt,[0.1,0.1],1)
divand.Quadtrees.add!(qt,[0.2,0.2],2)
divand.Quadtrees.add!(qt,[0.7,0.7],3)
divand.Quadtrees.add!(qt,[0.9,0.1],4)

divand.Quadtrees.split!(qt)


@test divand.Quadtrees.isleaf(qt) == false
@test divand.Quadtrees.isleaf(qt.children[1]) == true

X = rand(2,10000)
attribs = collect(1:size(X,2))
qt2 = divand.Quadtrees.QT(X,attribs)
divand.Quadtrees.rsplit!(qt2,5)

#rplot(qt2)
#plot(X[1,:],X[2,:],"b.")

xmin = [0.3,0.3]
xmax = [0.51,0.51]

function simplesearch(X,attribs,xmin,xmax)
    sel = trues(size(X,2))
    for j = 1:size(X,1)
        sel = sel .& (xmin[j] .<= X[j,:] .<= xmax[j])
    end
    ind = find(sel)
    return X[:,ind],attribs[ind]
end


xref,attribs_ref = simplesearch(X,attribs,xmin,xmax)

attribs_res = divand.Quadtrees.within(qt2,xmin,xmax)

#@show xref

@test sort(attribs_ref) == sort(attribs_res)


# progressively add all points
qt3 = divand.Quadtrees.QT(Int,[0.,0.],[1.,1.])

for i = 1:size(X,2)
    divand.Quadtrees.add!(qt3,X[:,i],i)
end

@test divand.Quadtrees.count(qt3) == size(X,2)

attribs_res = divand.Quadtrees.within(qt3,xmin,xmax)
@test sort(attribs_ref) == sort(attribs_res)

# Test in 1D - 4D

for n = 1:4
    #@show n
    X = rand(n,100)
    attribs = collect(1:size(X,2))

    qtND = divand.Quadtrees.QT(X,attribs)
    #@show qtND.points[1:5,:]

    divand.Quadtrees.rsplit!(qtND)

    @test divand.Quadtrees.ndims(qtND) == n
    @test divand.Quadtrees.count(qtND) == size(X,2)

    xmin = fill(0.0,(n,))
    xmax = fill(0.5,(n,))

    xref,attribs_ref = simplesearch(X,attribs,xmin,xmax)
    attribs_res = divand.Quadtrees.within(qtND,xmin,xmax)

    @test sort(attribs_ref) == sort(attribs_res)

    s = IOBuffer()
    show(s,qtND)
    @test contains(String(take!(s)),"Node")


end

# duplicates
x = Float64[i for i = 1:10, j = 1:11]
y = Float64[j for i = 1:10, j = 1:11]

value = x+y;
sel = [1:length(x);1];
dup = divand.Quadtrees.checkduplicates((x[sel],y[sel]),value[sel],[0.1,0.1],0.01)
@test dup == [[1,length(sel)]]

dup = divand.Quadtrees.checkduplicates((x[:],y[:]),value[:],
                                       (x[1:2],y[1:2]),value[1:2],[0.1,0.1],0.01)
@test dup == [[1],[2]]
