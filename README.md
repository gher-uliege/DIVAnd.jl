# DIVAnd
<div align="center"> <img src="docs/src/assets/logo.png"></img></div>

---

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Build Status](https://github.com/gher-uliege/DIVAnd.jl/workflows/CI/badge.svg)](https://github.com/gher-uliege/DIVAnd.jl/actions)
[![codecov](https://codecov.io/gh/gher-uliege/DIVAnd.jl/graph/badge.svg?token=iIDy1RvGXU)](https://codecov.io/gh/gher-uliege/DIVAnd.jl)
[![documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gher-uliege.github.io/DIVAnd.jl/stable/)
[![documentation latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://gher-uliege.github.io/DIVAnd.jl/latest/)
[![DOI](https://zenodo.org/badge/79277337.svg)](https://zenodo.org/badge/latestdoi/79277337)
![GitHub top language](https://img.shields.io/github/languages/top/gher-uliege/DIVAnd.jl)
<a href="https://archive.softwareheritage.org/swh:1:dir:d1dc29fdd3ca9b817a26a387dad5fc6f06060d37;origin=https://github.com/gher-uliege/DIVAnd.jl;visit=swh:1:snp:18bc3724a38165f88e21358f2f4b994cb2d548db;anchor=swh:1:rev:c071bebf0aafc211aaa3ef2ac0babd0ce80fc232">
<img src="https://archive.softwareheritage.org/badge/swh:1:dir:d1dc29fdd3ca9b817a26a387dad5fc6f06060d37/" alt="Archived | swh:1:dir:d1dc29fdd3ca9b817a26a387dad5fc6f06060d37"/>
</a>


`DIVAnd` (Data-Interpolating Variational Analysis in n dimensions) performs an n-dimensional variational analysis/gridding of arbitrarily located observations. Observations will be interpolated/analyzed on a curvilinear grid in 1, 2, 3 or more dimensions. In this sense it is a generalization of the original two-dimensional DIVA version (still available here https://github.com/gher-uliege/DIVA but not further developed anymore).

The method bears some similarities and equivalences with Optimal Interpolation or Krigging in that it allows to create a smooth and continous field from a collection of observations, observations which can be affected by errors. The analysis method is however different in practise, allowing to take into account topological features, physical constraints etc in a natural way. The method was initially developped with ocean data in mind, but it can be applied to any field where localized observations have to be used to produce gridded fields which are "smooth".

See also https://gher-uliege.github.io/DIVAnd-presentation/#1

Please cite this paper as follows if you use `DIVAnd` in a publication:

Barth, A., Beckers, J.-M., Troupin, C., Alvera-Azcárate, A., and Vandenbulcke, L.: DIVAnd-1.0: n-dimensional variational data analysis for ocean observations, Geosci. Model Dev., 7, 225-241, doi:[10.5194/gmd-7-225-2014](https://doi.org/10.5194/gmd-7-225-2014), 2014.

(click [here](./data/DIVAnd.bib) for the BibTeX entry).

# Summary of features

* N-Dimensional analysis/interpolation
* Scattered data
* Noise allowed
* Physical constraints can be added
* Inequality constraints can be added
* Topological constraints are handled naturally (barriers, holes)
* Analysis error maps can be estimated
* Periodicity in selected directions can be enforced
* Multivariate data can be used (experimental)
* The output grid can be curvilinear
* Instead of interpolating scattered data you can also peform Kernel Density Estimations with the points.


# Installing

You need [Julia](http://julialang.org) (version 1.6 or later) to run `DIVAnd`. The command line version is sufficient for `DIVAnd`.
Inside a Julia terminal, you can download and install the package by issuing:

```julia
using Pkg
Pkg.add("DIVAnd")
```

It is not recommended to download the source of `DIVAnd.jl` directly (using the green *Clone or Download* button above) because this by-passes Julia's package manager and you would need to install the dependencies of `DIVAnd.jl` manually.

## Cloud environement

DIVAnd is also available in the BlueCloud virtual research environement implemented by D4Science:
https://blue-cloud.d4science.org/group/coastalcurrentsfromobservations

Note that BlueCloud supports several authentication mechanisms, it is quite likely that you have already credentials that can be used to sign-in.
In `Analytics`, choose `JupyterLab on D4Science` then select the DIVAnd environement.

# Updating DIVAnd

To update DIVAnd, run the following command and restart Julia (or restart the jupyter notebook kernel using `Kernel` -> `Restart`):

```julia
using Pkg
Pkg.update("DIVAnd")
```

Note that Julia does not directly delete the previous installed version.
To check if you have the latest version run the following command:

```julia
using Pkg
Pkg.status()
```

The latest version number is available from [here](https://github.com/gher-uliege/DIVAnd.jl/releases).

To explicitly install a given version `X.Y.Z` you can also use:

```julia
using Pkg
Pkg.add(name="DIVAnd", version="X.Y.Z")
```
Or the master version:

```julia
using Pkg
Pkg.add(name="DIVAnd", rev="master")
```

# Testing

A test script is included to verify the correct functioning of the toolbox.
The script should be run in a Julia session.
Make sure to be in a directory with write-access (for example your home directory).
You can change the directory to your home directory with the `cd(homedir())` command.

```julia
using Pkg
Pkg.test("DIVAnd")
```

All tests should pass without error (it can take several minutes).

```
INFO: Testing DIVAnd
Test Summary: | Pass  Total
  DIVAnd      |  461    461
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

see also https://gher-uliege.github.io/DIVAnd.jl/latest/index.html

## Example

[DIVAnd_simple_example_4D.jl](https://github.com/gher-uliege/DIVAnd.jl/blob/master/examples/DIVAnd_simple_example_4D.jl) is a basic example in fours dimensions. The call to `DIVAndrun` looks like this:

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

More examples are available in the notebooks from the [Diva Workshop](https://github.com/gher-uliege/Diva-Workshops).

## Note on which analysis function to use

`DIVAndrun` is the core analysis function in n dimensions. It does not know anything about the physical parameters or units you work with. Coordinates can also be very general. The only constraint is that the metrics `(pm,pn,po,...)` when multiplied by the corresponding length scales `len` lead to non-dimensional parameters. Furthermore the coordinates of the output grid `(xi,yi,zi,...)` need to have the same units as the observation coordinates `(x,y,z,...)`.

`DIVAndfun` is a version with a minimal set of parameters (the coordinates and values of observations, i.e.  `(x,f)`, the remaining parameters being optional) and provides an interpolation *function* rather than an already gridded field. 

`diva3D` is a higher-level function specifically designed for climatological analysis of data on Earth, using longitude/latitude/depth/time coordinates and correlations length in meters. It makes the necessary preparation of metrics, parameter optimizations etc you normally would program yourself before calling the analysis function `DIVAndrun`.

`DIVAnd_heatmap` can be used for additive data and produces Kernel Density Estimations.


`DIVAndgo` is only needed for very large problems when a call to `DIVAndrun` leads to memory or CPU time problems. This function tries to decide which solver (direct or iterative) to use and how to make an automatic domain decomposition. Not all options from `DIVAndrun` are available.

If you want to try out multivariate approaches, you can look at `DIVAnd_multivarEOF` and `DIVAnd_multivarJAC`

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
Tools to help you are included in  ([DIVAnd_cv.jl](https://github.com/gher-uliege/DIVAnd.jl/blob/master/src/DIVAnd_cv.jl)).


## Note about the error fields

`DIVAnd` allows the calculation of the analysis error variance, scaled by the background error variance. Though it can be calculated "exactly" using the diagonal of the error covariance matrix s.P, it is generally too costly and approximations are provided. All of them are accessible as options via `DIVAnd_errormap` or you can let `DIVAnd` decide which version to use (possibly by specifying if you just need a quick estimate or a version closer the theoretical estimate) (see [Beckers et al 2014](https://doi.org/10.1175/JTECH-D-13-00130.1) )

## Advanced usage

### Additional constraint

An arbitrary number of additional quadratic constraints can be included to the cost function which should have the following form:

*J*(**x**) = ∑<sub>*i*</sub> (**C**<sub>*i*</sub> **x**  - **z**<sub>*i*</sub>)ᵀ **Q**<sub>*i*</sub><sup>-1</sup> (**C**<sub>*i*</sub> **x** - **z**<sub>*i*</sub>)

For every constrain, a structure with the following fields is passed to `DIVAnd`:

* `yo`: the vector **z**<sub>*i*</sub>
* `H`: the matrix **C**<sub>*i*</sub>
* `R`: the matrix **Q**<sub>*i*</sub> (symmetric and positive defined)

Internally the observations are also implemented as constraint defined in this way.

### Additional inequality constraint

An arbitrary number of additional inequality constraints can be included and which should have the following form:

(**H**<sub>*i*</sub> **x**  > **yo**<sub>*i*</sub>)

For every constraint, a structure with the following fields is passed to `DIVAnd`:

* `yo`: a vector
* `H`: a matrix



## Run notebooks on a server which has no graphical interface

On the server, launch the notebook with:
```bash
jupyter-notebook --no-browser --ip='0.0.0.0' --port=8888
```
where the path to `jupyter-notebook` might have to be adapted, depending on your installation. The `ip` and `port` parameters can also be modified.

Then from the local machine it is possible to connect to the server through the browser.

Thanks to Lennert and Bart (VLIZ) for this trick.

# Example data

Some examples in `DIVAnd.jl` use a quite large data set which cannot be efficiently distributed through `git`. This data can be downloaded from the URL https://dox.ulg.ac.be/index.php/s/Bo01EicxnMgP9E3/download. The zip file should be decompressed and the directory `DIVAnd-example-data` should be placed on the same level than the directory `DIVAnd.jl`.

# Reporting issues

Please include the following information when reporting an issue:

* Version of Julia
* Version of DIVAnd
* Operating system
* Full screen output preferably obtained by setting `ENV["JULIA_DEBUG"] = "DIVAnd"`.
* Full stack strace with error message
* A short description of the problem
* The command and their arguments which produced the error

Note that only [official julia builds](https://julialang.org/downloads/) are supported. 

In all cases, if we provide a tentative solution, please provide a feedback in all cases (whether it solved your issue or not).

# Fun

An [educational web application](http://data-assimilation.net/Tools/divand_demo/html/) has been developed to reconstruct a field based on point "observations". The user must choose in an optimal way the location of 10 observations such that the analysed field obtained by `DIVAnd` based on these observations is as close as possible to the original field.

# You do not want to use Julia

You should really reconsider and try out Julia. It is easy to use and provides the native interface to `DIVAnd`.

If you have a stable workflow using python, into which you want to integrate `DIVAnd`, you might try

https://github.com/gher-uliege/DIVAnd.py
