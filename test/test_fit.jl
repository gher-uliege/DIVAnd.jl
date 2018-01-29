using Base.Test
import divand

# test data for basic statistics
x = [1.,2.,3.,4.]
y = -[1.,2.,3.,4.]
sumx = sum(x)
sumx2 = sum(x.^2)

sumy = sum(y)
sumy2 = sum(y.^2)

sumxy = sum(x .* y)

meanx,stdx = divand.stats(sumx,sumx2,length(x))
meanx,meany,stdx,stdy,covar,corr = divand.stats(sumx,sumx2,sumy,sumy2,sumxy,length(x))

@test meanx ≈ mean(x)
@test stdx ≈ std(x)

@test corr ≈ -1

# DIVA fit

min_count = 50

# test data
nobs = 1000;
x = (2*pi*rand(nobs), 2*pi*rand(nobs))
v = sin.(x[1]) .* sin.(x[2]);

#@code_warntype fit(x,v,distbin,min_count)

# isotropic fit
len,var0,distx,covar,fitcovar = divand.fit_isotropic(x,v,distbin,min_count)
@test len < 2

# anisotropic fit
var0opt,lensopt,distx,covar,fitcovar = divand.fit(x,v,distbin,min_count)
@test all(lensopt .<= 2)
