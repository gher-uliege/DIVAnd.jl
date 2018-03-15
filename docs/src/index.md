
divand


# divand.jl documentation



```@docs
divand.diva3d
divand.divandrun
divand.divandgo
divand.divand_averaged_bg
divand.SDNMetadata
divand.save
divand.loadbigfile
divand.checkobs
divand.smoothfilter
divand.Anam.loglin
divand.Anam.logit
divand.divadoxml
divand.random
divand.distance
divand.interp
divand.backgroundfile
divand.Quadtrees.checkduplicates
```

# Bathymetry and spatial-temporal domain

```@docs
divand.load_bath
divand.extract_bath
divand.load_mask
divand.divand_metric
divand.domain
divand.divand_rectdom
divand.divand_squaredom
divand.TimeSelectorYW
divand.TimeSelectorYearListMonthList
```

# Load observations

```@docs
divand.saveobs
divand.loadobs
divand.NCSDN.load
divand.NCSDN.loadvar
divand.ODVspreadsheet.loaddata
divand.ODVspreadsheet.parsejd
```

# Parameter optimization

```@docs
divand.fit_isotropic
divand.fit
divand.divand_cv
divand.empiriccovar
divand.fithorzlen
divand.fitvertlen
divand.lengraddepth
divand.divand_cvestimator
divand.weight_RtimesOne
divand.Rtimesx!
```

# Vocabulary


```@docs
divand.Vocab.@urn_str
divand.Vocab.CFVocab
Base.haskey(collection::divand.Vocab.CFVocab,stdname)
divand.Vocab.SDNCollection
divand.Vocab.prefLabel
divand.Vocab.altLabel
divand.Vocab.notation
divand.Vocab.definition
divand.Vocab.resolve
divand.Vocab.find(c::divand.Vocab.Concept,name,collection)
divand.Vocab.description
divand.Vocab.canonical_units
divand.Vocab.splitURL
```

# Internal API or advanced usage


## State vector

```@docs
divand.statevector
divand.pack
divand.unpack
Base.sub2ind
Base.ind2sub
Base.length
```


## ODV files

```@docs
divand.ODVspreadsheet.listSDNparams
divand.ODVspreadsheet.load
divand.ODVspreadsheet.localnames
divand.ODVspreadsheet.Spreadsheet
divand.ODVspreadsheet.loadprofile
divand.ODVspreadsheet.loaddataqv
divand.ODVspreadsheet.SDNparse!
divand.ODVspreadsheet.colnumber
divand.ODVspreadsheet.nprofiles
```

## Operators

```@docs
divand.sparse_interp
divand.sparse_interp_g
divand.sparse_gradient
divand.sparse_diff
divand.matfun_trim
divand.matfun_stagger
divand.matfun_diff
divand.matfun_shift
```

## Quadtree

```@docs
divand.Quadtrees.QT
divand.Quadtrees.rsplit!
divand.Quadtrees.add!
divand.Quadtrees.within
divand.Quadtrees.bitget
divand.Quadtrees.inside
divand.Quadtrees.intersect
divand.Quadtrees.split!
```


## Conjugate gradient

```@docs
divand.conjugategradient
divand.pc_none!
divand.checksym
```

## Utility functions

```@docs
divand.divand_laplacian
divand.divand_obscovar
divand.divand_adaptedeps2
divand.divand_diagHKobs
divand.divand_residual
divand.divand_erroratdatapoints
divand.divand_iBpHtiRHx!
divand.divand_GCVKii
divand.divand_fittocpu
divand.divand_background
divand.divand_obs
divand.divand_bc_stretch
divand.divand_diagHK
divand.divand_kernel
divand.divand_residualobs
divand.divand_aexerr
divand.divand_cpme
divand.divand_Lpmnrange
divand.divand_pc_sqrtiB
divand.divand_pc_none
divand.divand_GCVKiiobs
divand.divand_cutter
divand.divand_qc
divand.divand_solve!
divand.divand_sampler
divand.divandjog
divand.divand_background_components
divand.stats
divand.statpos
divand.blkdiag
Base.findfirst
divand.formatsize
divand.interp!
divand.ufill
divand.jmBix
divand.cgradient
divand.fzero
divand.localize_separable_grid
divand.decompB!
divand.varanalysis
divand.len_harmonize
divand.alpha_default
divand.ncfile
divand.writeslice
divand.encodeWMSStyle
divand.loadoriginators
```



# Examples

To run the example, you need to install `PyPlot`.
In the folder `examples` of divand, you can run e.g. the example `divand_simple_example_1D.jl` by issuing:

```julia
# cd("/path/to/divand/examples")
include("divand_simple_example_1D.jl")
```

Replace `/path/to/divand/` by the installation directory of divand which is the output of `Pkg.dir("divand")` if you installed `divand` using Julias package manager.



# Information for developers

## Update the documentation

Install

```julia
Pkg.add("Documenter")
```

# Troubleshooting

## No plot windows

If the following command doesn't produce any figure
```julia
using PyPlot
plot(0, 1)
```
a possible solution is to modify the *backend*: this is done by editing the python configuration file
[matplotlibrc](http://matplotlib.org/users/customizing.html#the-matplotlibrc-file). The location of this file is obtained in python with:

```python
import matplotlib
matplotlib.matplotlib_fname
```
which, in my case, returns
```'~/.config/matplotlib/matplotlibrc'```
