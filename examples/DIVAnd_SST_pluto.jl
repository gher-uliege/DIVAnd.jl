### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 1c056adc-2e9a-11eb-3d81-0d29e98a3551
begin
    # We set up a new environment for this notebook
    import Pkg
    Pkg.activate(mktempdir())
    Pkg.add("PyCall")
    using PyCall
    try
        pyimport("matplotlib")
    catch
        run(`$(PyCall.python) -m pip install matplotlib`)
    end
    Pkg.add("PlutoUI")
    Pkg.add("DIVAnd")
    Pkg.add("PyPlot")
    using PlutoUI
    using DIVAnd
    using PyPlot
    using NCDatasets
	using Statistics
	using Random
end


# ╔═╡ 3e8cc514-2e9a-11eb-35ad-8f110251cb48
begin
     url = "http://psl.noaa.gov/thredds/dodsC/Datasets/noaa.oisst.v2/sst.ltm.1961-1990.nc"

    ds = NCDataset(url)
    vfull = reverse(nomissing(ds["sst"][:,:,1],NaN),dims=2)
    lon = ds["lon"][:]
    lat = reverse(ds["lat"][:])
    close(ds)

	mask = isfinite.(vfull)

	# remove zonal average
	v0 = copy(vfull)
	v0[.!mask] .= 0	
	v = vfull .- sum(v0,dims=1) ./ sum(mask,dims=1)
	
    sz  = size(v)
    xi,yi =  DIVAnd.ndgrid(lon,lat)

    i,j =  DIVAnd.ndgrid(1:sz[1],1:sz[2])
	
    # scale factor; inverse of the resolution
    pm = ones(sz) ./ ((xi[2,1]-xi[1,1]) .* cosd.(yi));
    pn = ones(sz) / (yi[1,2]-yi[1,1]);
	
	md"""### Illustration of DIVAnd

	We use the Reynolds et al. 2002 [OI SST](https://www.psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html) for the month January and remove the zonal average. We extract pseudo-observations at random location and aim to reconstruct the field from these data points.
	"""
end

# ╔═╡ b39c60e8-2e9b-11eb-3c60-7fe7082d59e2
@bind nobs  Slider(1:10000,default=1000)

# ╔═╡ b39c60e8-2e9b-11eb-3c60-7fe7082d59e3
md"Number of observations (`nobs` = $nobs)"

# ╔═╡ b39c60e8-2e9b-11eb-3c60-7fe7082d59e6
@bind noise  Slider(10. .^ (-2:0.1:2),default=0.1)

# ╔═╡ b39c60e8-2e9b-11eb-3c60-7fe7082d59e5
md"Noise on observation (`noise` = $noise)"

# ╔═╡ b39c60e8-2e9b-11eb-3c60-7fe7082d59e0
@bind len  Slider(LinRange(1,20,101),default=4)

# ╔═╡ b39c60e8-2e9b-11eb-3c60-7fe7082d59e1
md"Correlation length in arc degrees (1° ~ 111 km) (`len` = $len)"

# ╔═╡ 4482ce12-2e9c-11eb-028d-2100248f0f9a
@bind epsilon2  Slider(10. .^ (-2:0.1:2),default=0.1)

# ╔═╡ 2c3feeac-2f00-11eb-24bb-2b7015c0551e
md"Uncertainty of the observation (`epsilon2` = $epsilon2)"

# ╔═╡ b39c60e8-2e9b-11eb-3c60-7fe7082d59e4
@bind seed Slider(0:100,default=42)

# ╔═╡ eae323c7-b014-418d-80ad-bb752e67d5ac
begin
    Random.seed!(seed)
    indexobs = shuffle(findall(mask[:] .== 1))[1:nobs]
    vobs = v[indexobs] + noise * randn(nobs)
    vi,s = DIVAndrun(
        mask,(pm,pn),(xi,yi),(xi[indexobs],yi[indexobs]),
        vobs,
        len,epsilon2
    )
	close("all")
    figure(figsize = (7.8,5.6))
    function map(; cl = extrema(filter(isfinite,v)))
	axis("equal")
    xlim(extrema(lon))
    ylim(extrema(lat))
	clim(cl)
	colorbar(orientation="horizontal")
    contourf(xi,yi,mask,levels = [0,.5],cmap = "gray")
    end
	subplot(2,2,1)
    pcolormesh(xi,yi,v); title("True field")
	map()
	
	subplot(2,2,2)
	scatter(xi[indexobs],yi[indexobs],2,vobs,edgecolor="w",linewidth=0.);
	title("Observation")
    map()
	
	subplot(2,2,3)
    pcolormesh(xi,yi,vi); title("Analysis field")
	map()

	subplot(2,2,4)
    pcolormesh(xi,yi,vi - v); title("Analysis - true field")
	map(cl = (-1,1))
    
    gcf()
end

# ╔═╡ 4b9ddd5e-fd8f-4124-a1c2-3958fd926104
md"RMS difference between analysis and true field = $(sqrt(mean(filter(isfinite,(vi - v).^2))))"

# ╔═╡ b39c60e8-2e9b-11eb-3c60-7fe7082d59e7
md"Random seed (`seed` = $seed)"

# ╔═╡ 328a430d-1e7e-4ed9-8393-edb4e505ae0e


# ╔═╡ Cell order:
# ╟─1c056adc-2e9a-11eb-3d81-0d29e98a3551
# ╟─3e8cc514-2e9a-11eb-35ad-8f110251cb48
# ╟─b39c60e8-2e9b-11eb-3c60-7fe7082d59e3
# ╟─b39c60e8-2e9b-11eb-3c60-7fe7082d59e2
# ╟─b39c60e8-2e9b-11eb-3c60-7fe7082d59e5
# ╟─b39c60e8-2e9b-11eb-3c60-7fe7082d59e6
# ╟─b39c60e8-2e9b-11eb-3c60-7fe7082d59e1
# ╟─b39c60e8-2e9b-11eb-3c60-7fe7082d59e0
# ╟─2c3feeac-2f00-11eb-24bb-2b7015c0551e
# ╟─4482ce12-2e9c-11eb-028d-2100248f0f9a
# ╟─eae323c7-b014-418d-80ad-bb752e67d5ac
# ╟─4b9ddd5e-fd8f-4124-a1c2-3958fd926104
# ╟─b39c60e8-2e9b-11eb-3c60-7fe7082d59e7
# ╟─b39c60e8-2e9b-11eb-3c60-7fe7082d59e4
# ╟─328a430d-1e7e-4ed9-8393-edb4e505ae0e
