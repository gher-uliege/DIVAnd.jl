# Perform a divand benchmark
# A 2d-variational analysis is performed over a square domain
# as described in http://www.geosci-model-dev.net/7/225/2014/gmd-7-225-2014.pdf
# with different number of grid points (100 x 100, 200 x 200, ... 600 x 600).
# For every domain size, the benchmark is run 10 times and the median time is
# computed.
#
# Input (optional):
#   name: if name is present then the results will be saved in the file called
#     name
#
# Output:
#   median_time: the median of the run-time in seconds
#   ng: number of grid points along one dimension
#   time: 2-d array with all results
#
# To run this benchmark you need to install divand, for example
#
# pkg install -forge divand
# pkg load divand

# Alexander Barth
# GPLv2 or later

function test_2dvar_benchmark(name)

  @printf("Running divand benchmark in 2 dimensions\n");

  # domain sizes
  ng = 100:100:500;

  time = zeros(10,length(ng))
  RMS = zeros(10,length(ng))

  for j=1:10
    for i=1:length(ng)
      time[i,j],RMS[i,j] = benchmark2d(ng[i]);
      if (RMS[i,j] > 0.2)
        error("unexpected large RMS. Results might be wrong");
      end

      @printf("size %5d time %10.4f \n",ng[i],time[i,j]);
    end
  end

  @printf("\nMedian results\n");

  median_time = median(time,2);

  for i=1:length(ng)
    @printf("size %5d time %10.4f \n",ng[i],median_time[i]);
  end

  fname = "test_2dvar_benchmark_$(name).mat"
  @printf("save result in file %s\n",fname)
  matwrite(fname,Dict("time" => time,"RMS" => RMS,"ng" => ng))

  return median_time,ng,time
end

function benchmark2d_repeat(ng,ntimes)
    times = zeros(ntimes)
    RMS = zeros(ntimes)

    times[1],RMS[1] = benchmark2d(ng)

    for i = 1:ntimes
        times[i],RMS[1] = benchmark2d(ng)
    end

    mad(x) = median(abs(x - median(x)))

    stat = Dict{String,Any}([(string(f),f(times)) for f in [mean,std,median,mad,minimum,maximum]])
    stat["samples"] = length(times)
    stat["times"] = times
    stat["RMS"] = RMS
    return stat
end


function benchmark2d(ng)

  mg = ng/5;
  len = 10/ng;

  epsilon2 = 0.05;

  f(x,y) = cos(2*pi*ng*x/20) .* cos(2*pi*ng*y/20);


  # grid of background field
  xi,yi = ndgrid(linspace(0,1,ng),linspace(0,1,ng));
  vi = f(xi,yi);

  mask = trues(size(xi));
  pm = ones(size(xi)) / (xi[2,1]-xi[1,1]);
  pn = ones(size(xi)) / (yi[1,2]-yi[1,1]);

  # grid of observations
  x,y = ndgrid(linspace(1e-6,1-1e-6,mg),linspace(1e-6,1-1e-6,mg));
  x = x[:];
  y = y[:];
  v = f(x,y);

  t1 = time_ns()
  va,s = divandrun(mask,(pm,pn),(xi,yi),(x,y),v,len,epsilon2);
  t2 = time_ns()
  time = (t2 - t1)/1e9
  RMS = rms(va,vi);

  return time,RMS
end


function rms(x,y)

  d = x-y;

  m = !isnan.(d);
  r = mean(d[m].^2);

   r = sqrt.(r);
  return r
end
