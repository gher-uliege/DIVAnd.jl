# Testing the kernel of divand

using Base.Test

"""
mask,xyi,pmn = divand_squaredom(n,coord)

Create a "square" domain in `n` dimensions with the coordinates `coord`
assuming a Catersian metric. This functions returns
the mask `mask`, the coordinates `(xi,yi,...)` and the metric `(pm,pn...)`.

# Example

mask,(pm,pn),(xi,yi) = divand_squaredom(2,linspace(0,1,50))
"""
function divand_squaredom(n,coord)
    coords = ([coord for i = 1:n]...)
    return divand_rectdom(coords...)
end


"""
mask,xyi,pmn = divand_squaredom(n,coord)

Create a "square" domain in `n` dimensions with the coordinates `coord`
assuming a Catersian metric. This functions returns
the mask `mask`, the coordinates `(xi,yi,...)` and the metric `(pm,pn...)`.

# Example

mask,(pm,pn),(xi,yi) = divand_rectdom(linspace(0,1,50),linspace(0,1,50))
"""
function divand_rectdom(coords...)
    # grid of background field
    xyi = ndgrid(coords...)

    # mask (all points are valid)
    mask = trues(xyi[1])

    # metric (inverse of the resolution)
    pmn = ([ones(size(mask)) / (coords[i][2]-coords[i][1]) for i = 1:length(coords)]...)

    return mask,pmn,xyi
end


# correlation length
len = 0.2;

# value of the observations
f = [1.]

# normalized error variance
epsilon2 = 1.;

# dimension
for n = 1:3
    # domain
    mask,pmn,xyi = divand_squaredom(n,linspace(0,1,20))

    # grid of observations
    xy = ([[0.5] for i = 1:n]...)

    # make the analysis
    fi,s = divandrun(mask,pmn,xyi,xy,f,len,epsilon2);

    #@show maximum(fi),n,ndims(mask)

    if n < 3
        @test 0.48 <= maximum(fi) <= 0.52
    else
        @test 0.4 <= maximum(fi) <= 0.6
    end
end

# additional test in two dimensions
n = 2

# domain
mask,pmn,xyi = divand_squaredom(n,linspace(0,1,50))

# grid of observations
xy = ([[0.5] for i = 1:n]...)

# make the analysis
fi,s = divandrun(mask,pmn,xyi,xy,f,len,epsilon2,alpha=[2,4,2]);
@test 0.48 <= maximum(fi) <= 0.52

fi,s = divandrun(mask,pmn,xyi,xy,f,len,epsilon2,alpha=[1,0,1]);
@test 0.4 <= maximum(fi) <= 0.6
