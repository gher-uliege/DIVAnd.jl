# Testing the kernel of DIVAnd

using Base.Test


# correlation length
len = 0.2;

# value of the observations
f = [1.]

# normalized error variance
epsilon2 = 1.;

# dimension
for n = 1:3
    # domain
    mask,pmn,xyi = DIVAnd_squaredom(n,linspace(0,1,20))

    # grid of observations
    xy = ([[0.5] for i = 1:n]...)

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
mask,pmn,xyi = DIVAnd_squaredom(n,linspace(0,1,50))

# grid of observations
xy = ([[0.5] for i = 1:n]...)

# make the analysis
fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,alpha=[2,4,2]);
@test 0.47 <= maximum(fi) <= 0.53

fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,alpha=[1,0,1]);
@test 0.4 <= maximum(fi) <= 0.6
