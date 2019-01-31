function hmerge(f,L)
    # L ∼ (α * nmax)²
    # nmax ∼ √(L)/α

    weight0 = Float64.(isfinite.(f));

    mask,pmn = DIVAnd.DIVAnd_rectdom(1:size(f,1),1:size(f,2))
    ivol,nus = DIVAnd.DIVAnd_laplacian_prepare(mask,pmn,(ones(size(mask)),ones(size(mask))))

    α = 0.1;
    nmax = round(Int,sqrt(L)/α)
    @show nmax

    #nmax = 20;
    weight = similar(weight0);

    for k = 1:size(weight,3)
        wk0 = @view weight0[:,:,k]
        wk = @view weight[:,:,k]

        diffusionfix!(ivol,nus,α,nmax,wk0,wk)
    end
    f[.!isfinite.(f)] .= 0

    weight = weight.^2;
    #pcolor(copy(weight[:,:,1]'))
    
    f2 = (sum(weight .* f, dims = 3) ./ sum(weight, dims = 3))[:,:,1]
    return f2
end

function diffusionfix!(ivol,nus,α,nmax,x0,x)
    work1 = similar(x)
    x[:] = x0

    for niter = 1:nmax
        DIVAnd.DIVAnd_laplacian_apply!(ivol,nus,x,work1)
        for i in 1:length(x0)
           if x0[i] != 0
              x[i] = x[i] + α * work1[i]
           end
         end
    end

end

sz = (80,90,3)
f = fill(NaN,sz);
f[24:40,30:65,1] .= 1;
f[30:50,60:75,2] .= 2;
f[20:35,20:55,3] .= 1.5;

L = 4
f2 = hmerge(f,L);

#pcolor(copy(isfinite.(weight)'));
#pcolor(copy(isfinite.(weight0)'));

#pcolor(copy(weight[:,:,1]'));

pcolor(copy(f2'));
