# Testing the kernel of DIVAnd

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
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
    mask,pmn,xyi = DIVAnd_squaredom(n,Compat.range(0, stop = 1, length = 20))

    # grid of observations
    xy = ntuple(i -> [0.5], n)

    # make the analysis
    fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2);

    #@show maximum(fi),n,ndims(mask)

    if n < 3
        @test 0.45 <= maximum(fi) <= 0.55
    else
        @test 0.4 <= maximum(fi) <= 0.6
    end
end

# additional test in two dimensions
n = 2

# domain
mask,pmn,xyi = DIVAnd_squaredom(n,Compat.range(0, stop = 1, length = 50))

# grid of observations
xy = ntuple(i -> [0.5], n)

# make the analysis
fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,alpha=[2,4,2]);
@test 0.47 <= maximum(fi) <= 0.53

fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,alpha=[1,0,1]);
@test 0.4 <= maximum(fi) <= 0.6
