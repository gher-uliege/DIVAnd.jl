# Testing the kernel of DIVAnd

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
using DIVAnd

# correlation length
len = 0.2;

# value of the observations
f = [1.]

# normalized error variance
epsilon2 = 1.;

# dimension
for n = 1:3
    # domain
    mask,pmn,xyi = DIVAnd_squaredom(n,range(0, stop = 1, length = 20))

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
mask,pmn,xyi = DIVAnd_squaredom(n,range(0, stop = 1, length = 50))

# grid of observations
xy = ntuple(i -> [0.5], n)

# make the analysis
fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,alpha=[2,4,2]);
@test 0.47 <= maximum(fi) <= 0.53

fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,alpha=[1,0,1]);
@test 0.4 <= maximum(fi) <= 0.6



fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,primal=true,coeff_derivative2 = [1.,1.],coeff_laplacian = [0.,0.]);
@show maximum(fi)


# dimension
n = 3
mask,pmn,xyi = DIVAnd_squaredom(n,range(0, stop = 1, length = 30))

# grid of observations
xy = ntuple(i -> [0.5], n)

fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,primal=true);
@show maximum(fi)

# make the analysis
fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,primal=true,coeff_derivative2 = [1.,1.,1.],coeff_laplacian = [0.,0.,0.]);
@show maximum(fi)

fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,primal=true,coeff_derivative2 = [0.,0.,1.],coeff_laplacian = [1.,1.,0.]);
@show maximum(fi)


fi,s = DIVAndrun(mask,pmn,xyi,xy,f,len,epsilon2,primal=true,coeff_derivative2 = [0.,0.,0.1],coeff_laplacian = [1.,1.,0.9]);
@show maximum(fi)

#@test all(abs.(va[isfinite.(va)]) .< 1)
