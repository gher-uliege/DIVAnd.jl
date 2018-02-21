using Base.Test
import divand
using PyPlot


fname = joinpath(dirname(@__FILE__),"..","..","divand-example-data","BlackSea","Salinity.bigfile")


function vertcorrlen(x,value::Vector{T},z;
  distbin = collect([0.:50:400; 500:100:1000]),
  mincount = 50,
  maxpoints = 10000,
  nmean = 10,
  len0 = 50.,
  maxlen = 1000.,
  smoothz = 100.,
  ) where T

  function distfun(xi,xj)
    sqrt(sum(abs2,xi-xj))
  end

  zlevel2 = zero(T)

  function pickone(mask)
    ii = find(mask)
    jindex = rand(1:length(ii))
    return ii[jindex]
  end



  function vertchoose(x,zlevel)
    mask = abs.(zlevel - x[3]) .< 20
    j = pickone(mask);

    mask = falses(size(x[1]))
    for k = 1:length(x[1])
      mask[k] = distfun([x[1][k],x[2][k]],[x[1][j],x[2][j]]) < 0.5
    end

    jindex = pickone(mask)

    #@show j,jindex
    return j,jindex
  end


  function vchoose(dummy)
    vertchoose(x,zlevel2::Float64)
  end

  #@code_warntype vertchoose(x,zlevel2)
  @code_warntype vchoose(x)

  pmax = length(distbin)-1
  kmax = length(z)
  len = zeros(kmax)
  var0 = zeros(kmax)
  fitcovar = Array{T,2}(pmax,kmax)
  covar = Array{T,2}(pmax,kmax)
  distx = Vector{T}(pmax)


  # @time distx,covar,corr,varx,count = divand.empiriccovarmean(
  #      (depth[sel],),value[sel],distbin,mincount;
  #      maxpoints = maxpoints,
  #      nmean = nmean,
  #      choose = vchoose,
  #      distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj)))


  # k = 1
  # @time var0[k],len[k],distx,covar,fitcovar = divand.fit_isotropic(
  #      (depth[sel],),value[sel],distbin,mincount;
  #      maxpoints = maxpoints,
  #      nmean = nmean,
  #      maxlen = maxlen,
  #      choose = vchoose)


  # clf()
  #
  # k = 10
  # distx,covar,corr,varx,count = divand.empiriccovarmean(
  #     (depth[sel],),value[sel],distbin,mincount;
  #     maxpoints = maxpoints,
  #     nmean = nmean,
  #     choose = vchoose)
  #
  #     @time var0[k],len[k],distx,covar,fitcovar = divand.fit_isotropic(
  #         (depth[sel],),value[sel],distbin,mincount;
  #         maxpoints = maxpoints,
  #         nmean = nmean,
  #         maxlen = maxlen,
  #         len = len0,
  #         choose = vchoose)


  for k = 1:length(z)
    zlevel2 = Float64(z[k])
    @show zlevel2,maxpoints,nmean

    var0[k],len[k],distx[:],covar[:,k],fitcovar[:,k] = divand.fit_isotropic(
    (x[3],),value,distbin,mincount;
    maxpoints = maxpoints,
    nmean = nmean,
    maxlen = maxlen,
    len = len0,
    choose = vchoose)

     println("Data points at z=$(z[k]): $(length(v)), correlation length: $(len[k])")
    #plot(distx,covar, label = "empirical covariance")
    #plot(distx,fitcovar, label = "fitted function")
  end
  lenf = copy(len)
  if smoothz > 0
    divand.smoothfilter!(z,lenf,smoothz)
  end

  return lenf,Dict(
  :var0 => var0,
  :len => len,
  :distx => distx,
  :covar => covar,
  :fitcovar => fitcovar)

end



value,lon,lat,depth,time,ids = divand.loadbigfile(fname)
sel = Dates.month.(time) .== 1;
x = (lon[sel],lat[sel],depth[sel]);
v = value[sel]
z = [0.,10,100,200,300,400,500,700,1000,1500]
#z = [0.,10,100]


lenz,dbinfo = vertcorrlen(x,v,z)
