using Base.Test
import Quadtrees
    @testset "quadtree" begin

        srand(123)

        X = [0  0;
             1  0;
             1  1;
             0  1;
             0  0]

        qt = Quadtrees.QT(X,collect(1:size(X,1)))
        attribs_res = Quadtrees.within(qt,[0,0],[0.1,0.1])
        @test attribs_res == [1,5]


        @test [Quadtrees.bitget(42,i) for i = 6:-1:1] == [true,false,true,false,true,false]

        @test Quadtrees.inside([0,0],[1,1],[0.5,0.5]) == true
        @test Quadtrees.inside([0,0],[1,1],[1.5,1.5]) == false


        @test Quadtrees.intersect([0,0],[1,1],[0.5,0.5],[2,2]) == true
        @test Quadtrees.intersect([0,0],[1,1],[1.5,1.5],[2,2]) == false
        # one rectange contains the other
        @test Quadtrees.intersect([0,0],[1,1],[-1,-1],[2,2]) == true
        #@test_throws ArgumentError intersect([0,0],[1,1],[0.5,0.5],[2,2,3])

        qt = Quadtrees.QT(Int,[0.,0.],[1.,1.])
        Quadtrees.add!(qt,[0.1,0.1],1)
        Quadtrees.add!(qt,[0.2,0.2],2)
        Quadtrees.add!(qt,[0.7,0.7],3)
        Quadtrees.add!(qt,[0.9,0.1],4)

        Quadtrees.split!(qt)


        @test Quadtrees.isleaf(qt) == false
        @test Quadtrees.isleaf(qt.children[1]) == true

        X = rand(10000,2)
        attribs = collect(1:size(X,1))
        qt2 = Quadtrees.QT(X,attribs)
        Quadtrees.rsplit!(qt2,5)

        #rplot(qt2)
        #plot(X[:,1],X[:,2],"b.")

        xmin = [0.3,0.3]
        xmax = [0.51,0.51]

        function simplesearch(X,attribs,xmin,xmax)
            sel = trues(size(X,1))
            for j = 1:size(X,2)
                sel = sel .& (xmin[j] .<= X[:,j] .<= xmax[j])
            end
            ind = find(sel)
            return X[ind,:],attribs[ind]
        end


        xref,attribs_ref = simplesearch(X,attribs,xmin,xmax)

        attribs_res = Quadtrees.within(qt2,xmin,xmax)

        #@show xref

        @test sort(attribs_ref) == sort(attribs_res)


        # progressively add all points
        qt3 = Quadtrees.QT(Int,[0.,0.],[1.,1.])

        for i = 1:size(X,1)
            Quadtrees.add!(qt3,X[i,:],i)
        end

        @test Quadtrees.count(qt3) == size(X,1)

        attribs_res = Quadtrees.within(qt3,xmin,xmax)
        @test sort(attribs_ref) == sort(attribs_res)

        # Test in 1D - 4D

        for n = 1:4
            #@show n
            X = rand(100,n)
            attribs = collect(1:size(X,1))

            qtND = Quadtrees.QT(X,attribs)
            #@show qtND.points[1:5,:]

            Quadtrees.rsplit!(qtND)

            @test Quadtrees.ndims(qtND) == n
            @test Quadtrees.count(qtND) == size(X,1)

            xmin = fill(0.0,(n,))
            xmax = fill(0.5,(n,))

            xref,attribs_ref = simplesearch(X,attribs,xmin,xmax)
            attribs_res = Quadtrees.within(qtND,xmin,xmax)

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
        dup = Quadtrees.checkduplicates((x[sel],y[sel]),value[sel],[0.1,0.1],0.01)
        @test dup == [[1,length(sel)]]

        dup = Quadtrees.checkduplicates((x,y),value,(x[1:2],y[1:2]),value[1:2],[0.1,0.1],0.01)
        @test dup == [[1],[2]]

    end
