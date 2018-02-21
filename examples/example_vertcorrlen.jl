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
    zindex = Vector{Int}(length(value))
    nzindex = 0
        
  function pickone(mask)
    ii = find(mask)
    jindex = rand(1:length(ii))
    return ii[jindex]
  end



  function vertchoose(x,zlevel)
    #mask = abs.(zlevel - x[3]) .< 20
      #j = pickone(mask);
      maxntries = 10000
      maxntries2 = 1000
      j = -1
      jindex = -1
      
    for ntries = 1:maxntries  
      
    j = zindex[rand(1:length(zindex))] :: Int

    #mask = falses(size(x[1]))
    #for k = 1:length(x[1])
    #  mask[k] = distfun([x[1][k],x[2][k]],[x[1][j],x[2][j]]) < 0.5
      #end
    #jindex = pickone(mask)
      jindex = -1

        for ntries2 = 1:maxntries2
          k = rand(1:length(x[1]))
          if distfun([x[1][k],x[2][k]],[x[1][j],x[2][j]]) < 0.5
              jindex = k
              break
          end
        end

        if jindex != -1
            break
        end
    end

      if (j == -1) || (jindex == -1)
          error("fail to find enought pairs at z = $(zlevels2)")
      end
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
    zindex = find(abs.(zlevel2 - x[3]) .< 20)
      
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
z = [0.,10]


srand(12345);
@time lenz,dbinfo = vertcorrlen(x,v,z)

srand(12345);
@time lenz,dbinfo = vertcorrlen(x,v,z)


#=
srand(12345); @time lenz,dbinfo = vertcorrlen(x,v,z)

Data points at z=0.0: 3065, correlation length: 15.31864937639221
(zlevel2, maxpoints, nmean) = (10.0, 10000, 10)
Data points at z=10.0: 3065, correlation length: 127.22510065919151
184.761368 seconds (4.73 G allocations: 230.815 GiB, 25.09% gc time)
([71.2719, 71.2719], Dict{Symbol,Any}(Pair{Symbol,Any}(:fitcovar, [2.31656 0.126384; 0.0628743 0.105038; … ; 4.74498e-31 0.000139351; 6.87451e-35 5.24698e-5]),Pair{Symbol,Any}(:covar, [2.48498 1.59739; 0.0538359 0.047323; … ; 0.00147103 0.00147502; 0.00126514 0.00189619]),Pair{Symbol,Any}(:len, [15.3186, 127.225]),Pair{Symbol,Any}(:distx, [25.0, 75.0, 125.0, 175.0, 225.0, 275.0, 325.0, 375.0, 450.0, 550.0, 650.0, 750.0, 850.0, 950.0]),Pair{Symbol,Any}(:var0, [6.70391, 0.130344])))


=#
