# DIVAnd

[![Build Status Linux and macOS](https://travis-ci.org/gher-ulg/DIVAnd.jl.svg?branch=master)](https://travis-ci.org/gher-ulg/DIVAnd.jl)
[![Build Status Windows](https://ci.appveyor.com/api/projects/status/github/gher-ulg/DIVAnd.jl?branch=master&svg=true)](https://ci.appveyor.com/project/Alexander-Barth/DIVAnd-jl)

[![Coverage Status](https://coveralls.io/repos/gher-ulg/DIVAnd.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/gher-ulg/DIVAnd.jl?branch=master)
[![codecov.io](http://codecov.io/github/gher-ulg/DIVAnd.jl/coverage.svg?branch=master)](http://codecov.io/github/gher-ulg/DIVAnd.jl?branch=master)

<!--[![documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gher-ulg.github.io/DIVAnd.jl/stable/)-->
[![documentation latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://gher-ulg.github.io/DIVAnd.jl/latest/)

[![DOI](https://zenodo.org/badge/79277337.svg)](https://zenodo.org/badge/latestdoi/79277337)

`DIVAnd` (Data-Interpolating Variational Analysis in n dimensions) performs an n-dimensional variational analysis of arbitrarily located observations. Observations will be interpolated on a curvilinear grid in 2, 3 or more dimensions.

Please cite this paper as follows if you use `DIVAnd` in a publication:

Barth, A., Beckers, J.-M., Troupin, C., Alvera-Azcárate, A., and Vandenbulcke, L.: DIVAnd-1.0: n-dimensional variational data analysis for ocean observations, Geosci. Model Dev., 7, 225-241, doi:[10.5194/gmd-7-225-2014](http://dx.doi.org/10.5194/gmd-7-225-2014), 2014.

(click [here](./data/DIVAnd.bib) for the BibTeX entry).


# Installing

Under Linux you will also need the packages `make`, `gcc` and `netcdf` which you can install under Debian/Ubuntu with:

```bash
apt-get install make gcc libnetcdf-dev netcdf-bin
```

You need [Julia](http://julialang.org) (version 1.0 or 1.1) to run `DIVAnd`. The command line version is sufficient for `DIVAnd`.
Inside Julia, you can download and install the package by issuing:

```julia
using Pkg
Pkg.add(PackageSpec(name="DIVAnd", rev="master"))
```

For Julia 0.6, you can use the following:
```julia
Pkg.clone("https://github.com/gher-ulg/DIVAnd.jl") # only for Julia 0.6
```

It is not recommended to download the source of `DIVAnd.jl` directly (using the green *Clone or Download* button above) because this by-passes Julia's package manager and you would need to install the dependencies of `DIVAnd.jl` manually.


# Updating DIVAnd

To update DIVAnd, run the following command and restart Julia (or restart the jupyter notebook kernel):

```julia
Pkg.update()
```


# Testing

A test script is included to verify the correct functioning of the toolbox.
The script should be run in a Julia session.
Make sure to be in a directory with write-access (for example your home directory).
You can change the directory to your home directory with the `cd(homedir())` command.

```julia
Pkg.test("DIVAnd")
```

All tests should pass without error.

```
INFO: Testing DIVAnd
Test Summary: | Pass  Total
  DIVAnd      |   100     100
INFO: DIVAnd tests passed
```

The test suite will download some sample data.
You need to have Internet access and run the test function from a directory with write access.

# Documentation

The main routine of this toolbox is called `DIVAnd` which performs an n-dimensional variational analysis of arbitrarily located observations. Type the following in Julia to view a list of parameters:

```julia
using DIVAnd
?DIVAndrun
```

see also https://gher-ulg.github.io/DIVAnd.jl/latest/index.html

## Example

[DIVAnd_simple_example_4D.jl](https://github.com/gher-ulg/DIVAnd.jl/blob/master/examples/DIVAnd_simple_example_4D.jl) is a basic example in fours dimensions. The call to `DIVAndrun` looks like this:

```julia
fi,s = DIVAndrun(mask,(pm,pn,po,pq),(xi,yi,zi,ti),(x,y,z,t),f,len,epsilon2);

```
where
`mask` is the land-sea mask, usually obtained from the bathymetry/topography,
`(pm,pn,po,pq)` is a *n*-element tuple (4 in this case) containing the scale factors of the grid,
`(xi,yi,zi,ti)` is a *n*-element tuple containing the coordinates of the final grid,
`(x,y,z,t)` is a *n*-element tuple containing the coordinates of the observations,
`f` is the data anomalies (with respect to a background field),
`len` is the correlation length and
`epsilon2` is the error variance of the observations.

The call returns `fi`, the analyzed field on the grid `(xi,yi,zi,ti)`.

## Note on which analysis function to use

`DIVAndrun` is the core analysis function in n dimensions. It does not know anything about the physical parameters or units you work with. Coordinates can also be very general. The only constraint is that the metrics `(pm,pn,po,...)` when multiplied by the corresponding length scales `len` lead to non-dimensional parameters. Furthermore the coordinates of the output grid `(xi,yi,zi,...)` need to have the same units as the observation coordinates `(x,y,z,...)`.

`DIVAndgo` is only needed for very large problems when a call to `DIVAndrun` leads to memory or CPU time problems. This function tries to decide which solver (direct or iterative) to use and how to make an automatic domain decomposition. Not all options from `DIVAndrun` are available.

`diva3D` is a higher-level function specifically designed for climatological analysis of data on Earth, using longitude/latitude/depth/time coordinates and correlations length in meters. It makes the necessary preparation of metrics, parameter optimizations etc you normally would program yourself before calling the analysis function `DIVAndrun`.


## Note about the background field

If zero is not a valid first guess for your variable (as it is the case for e.g. ocean temperature), you have to subtract the first guess from the observations before calling `DIVAnd` and then add the first guess back in.

## Determining the analysis parameters

The parameter `epsilon2` and parameter `len` are crucial for the analysis.

`epsilon2` corresponds to the inverse of the [signal-to-noise ratio](https://en.wikipedia.org/wiki/Signal-to-noise_ratio). `epsilon2` is the normalized variance of observation error (i.e. divided by the background error variance). Therefore, its value depends on how accurate and how representative the observations are.
`len` corresponds to the correlation length and the value of `len` can sometimes be determined by physical arguments. Note that there should be one correlation length per dimension of the analysis.

One statistical way to determine the parameter(s) is to do a [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29).

1. choose, at random, a relatively small subset of observations (about 5%). This is the validation data set.
2. make the analysis without your validation data set
3. compare the analysis to your validation data set and compute the RMS difference
4. repeat steps 2 and 3 with different values of the parameters and try to minimize the RMS difference.

You can repeat all steps with a different validation data set to ensure that the optimal parameter values are robust.
Tools to help you are included in  ([DIVAnd_cv.jl](https://github.com/gher-ulg/DIVAnd.jl/blob/master/src/DIVAnd_cv.jl)).

## Advanced usage

### Additional constraint

An arbitrary number of additional constraints can be included to the cost function which should have the following form:

*J*(**x**) = ∑<sub>*i*</sub> (**C**<sub>*i*</sub> **x**  - **z**<sub>*i*</sub>)ᵀ **Q**<sub>*i*</sub><sup>-1</sup> (**C**<sub>*i*</sub> **x** - **z**<sub>*i*</sub>)

For every constrain, a structure with the following fields is passed to `DIVAnd`:

* `yo`: the vector **z**<sub>*i*</sub>
* `H`: the matrix **C**<sub>*i*</sub>
* `R`: the matrix **Q**<sub>*i*</sub> (symmetric and positive defined)

Internally the observations are also implemented as constraint defined in this way.

## Run notebooks on a server which has no graphical interface

On the server, launch the notebook with:
```bash
~/.julia/v0.6/Conda/deps/usr/bin/jupyter-notebook --no-browser --ip='0.0.0.0' --port=8888
```
where the path to `jupyter-notebook` might have to be adapted, depending on your installation. The `ip` and `port` parameters can also be modified.

Then from the local machine it is possible to connect to the server through the browser.

Thanks to Lennert and Bart (VLIZ) for this trick.

# Example data

Some examples in `DIVAnd.jl` use a quite large data set which cannot be efficiently distributed through `git`. This data can be downloaded from the URL https://dox.ulg.ac.be/index.php/s/Bo01EicxnMgP9E3/download. The zip file should be decompressed and the directory `DIVAnd-example-data` should be placed on the same level than the directory `DIVAnd.jl`.


# Fun

An [educational web application](http://data-assimilation.net/Tools/divand_demo/html/) has been developed to reconstruct a field based on point "observations". The user must choose in an optimal way the location of 10 observations such that the analysed field obtained by `DIVAnd` based on these observations is as close as possible to the original field.
