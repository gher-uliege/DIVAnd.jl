# divand



[![Build Status](https://travis-ci.org/gher-ulg/divand.jl.svg?branch=master)](https://travis-ci.org/gher-ulg/divand.jl)
[![Coverage Status](https://coveralls.io/repos/gher-ulg/divand.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/gher-ulg/divand.jl?branch=master) [![codecov.io](http://codecov.io/github/gher-ulg/divand.jl/coverage.svg?branch=master)](http://codecov.io/github/gher-ulg/divand.jl?branch=master)


`divand` performs an n-dimensional variational analysis of arbitrarily located observations. Observations will be interpolated on a curvilinear grid in 2, 3 or more dimensions.

Please cite this paper as follows if you use `divand` in a publication:

Barth, A., Beckers, J.-M., Troupin, C., Alvera-Azcárate, A., and Vandenbulcke, L.: divand-1.0: n-dimensional variational data analysis for ocean observations, Geosci. Model Dev., 7, 225-241, [doi:10.5194/gmd-7-225-2014](http://dx.doi.org/10.5194/gmd-7-225-2014), 2014.

# Installing

Your need [Julia](http://julialang.org) to run `divand`. The command line version is sufficient for `divand`.
Inside Julia, you can download and install the package by issuing:

```julia
Pkg.clone("https://github.com/gher-ulg/divand.jl")
```

# Testing

A test script is included to verify the correct functioning of the toolbox.
The script should be run in a Julia session.

```julia
Pkg.test("divand")
```

All tests should pass without warning or error.

```
INFO: Testing divand
Test Summary: | Pass  Total
  divand      |   37     37
INFO: divand tests passed
```


# Documentation

The main routine of this toolbox is called `divand` which performs an n-dimensional variational analysis of arbitrarily located observations. Type the following in Julia to view a list of parameters:

```julia
using divand
?divandrun
```

## Note

If zero is not a valid first guess for your variable (as it is the case for e.g. ocean temperature), you have to subtract the first guess from the observations before calling divand and then add the first guess back in.


## Advanced usage

### Additional constraint

An arbitrary number of additional constraints can be included to the cost function which should have the following form:


*J*(**x**) = ∑<sub>*i*</sub> (**C**<sub>*i*</sub> **x**  - **z**<sub>*i*</sub>)ᵀ **Q**<sub>*i*</sub><sup>-1</sup> (**C**<sub>*i*</sub> **x** - **z**<sub>*i*</sub>)

For every constrain, a structure with the following fields is passed to `divand`:

* `yo`: the vector **z**<sub>*i*</sub>
* `H`: the matrix **C**<sub>*i*</sub>
* `R`: the matrix **Q**<sub>*i*</sub> (symmetric and positive defined)

Internally the observations are also implemented as an additional constraint.

## Example

See the file [divand_simple_example.jl](https://github.com/gher-ulg/divand.jl/blob/master/examples/divand_simple_example.jl).


## Determining the parameters

The parameter `epsilon2` and parameter `len` are crucial for the analysis. `epsilon2` corresponds to the inverse of the [signal-to-noise ratio](https://en.wikipedia.org/wiki/Signal-to-noise_ratio). `epsilon2` is the normalizd variance of observation error (i.e. divided by the background error variance). Therefore, its value depends on how accurate and how representative the observations are. `len` corresponds to the correlation length and its value of `len` can sometimes be determined by physical arguments.
One statistical way to determine the parameter(s) is to do a [cross-validation](https://en.wikipedia.org/wiki/Cross-validation_%28statistics%29).

1. choose, at random, a relatively small subset of observations (about 5%). This is the validation data set.
2. make the analysis without your validation data set
3. compare the analysis to your validation data set and compute the RMS difference
4. repeat steps 2 and 3 with different values of the parameters and try to minimize the RMS difference.

You can repeat all steps with a different validation data set to ensure that the optimal parameter values are robust.
Tools to help you are included ([divand_cv.jl]).

# Fun

A [educational web application](http://data-assimilation.net/Tools/divand_demo/html/) has been developed to reconstruct a field based on point "observations". The user must choose in an optimal way the location of 10 observations such that the analysed field obtained by divand based on these observations is as close as possible to the original field.
