
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
divand.ODVspreadsheet.myparse
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
divand.divand_addc
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
divand.divand_cpme_go
divand.divand_datainboundingbox
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


If the installation of a package fails, it is recommended to update the local copy of the package list by issuing `Pkg.update()` to make sure that Julia knows about the latest version of these packages and then to re-try the installation of the problematic package. 
Julia calls the local copy of the packge list `METADATA`.
For example to retry the installation of EzXML issue the following command:

```julia
Pkg.update()
Pkg.add("EzXML")
```



## No plotting window appears

If the following command doesn't produce any figure
```julia
using PyPlot
plot(1:10)
```
A possible solution is to modify the *backend*: this is done by editing the python configuration file
[matplotlibrc](http://matplotlib.org/users/customizing.html#the-matplotlibrc-file). The location of this file is obtained in python with:

```python
import matplotlib
matplotlib.matplotlib_fname
```

Under Linux, this returns ```'~/.config/matplotlib/matplotlibrc'```.
To use the `TkAgg` backend, add the following to the file:

```
backend      : TkAgg
```

The `matplotlibrc` need to be created if it does not exists.

## Julia cannot connect to GitHub on Windows 7 and Windows Server 2012

Cloning METADATA or downloading a julia packages fails with:

```
GitError(Code:ECERTIFICATE, Class:OS, , user cancelled certificate checks: )
```

The problem is that Windows 7 and Windows Server 2012 uses outdated encryption protocols. The solution is to run the 
"Easy fix" tool from the [Microsoft support page](https://stackoverflow.com/questions/49065986/installation-of-julia-on-windows7-64-bit)

## MbedTLS.jl does not install on Windows 7


The installion of `MbedTLS.jl` fails with the error message:

```
INFO: Building MbedTLS                                                                                                                                    
Info: Downloading https://github.com/quinnj/MbedTLSBuilder/releases/download/v0.6/MbedTLS.x86_64-w64-mingw32.tar.gz to C:\Users\Jeremy\.julia\v0.6\MbedTLS
\deps\usr\downloads\MbedTLS.x86_64-w64-mingw32.tar.gz...                                                                                                  
Exception setting "SecurityProtocol": "Cannot convert null to type "System.Net.SecurityProtocolType" due to invalid enumeration values. Specify one of th 
e following enumeration values and try again. The possible enumeration values are "Ssl3, Tls"."                                                           
At line:1 char:35                                                                                                                                         
+ [System.Net.ServicePointManager]:: <<<< SecurityProtocol =                                                                                              
    + CategoryInfo          : InvalidOperation: (:) [], RuntimeException                                                                                  
    + FullyQualifiedErrorId : PropertyAssignmentException                                                                                                 
    [...]
```

See also the issue https://github.com/JuliaWeb/MbedTLS.jl/issues/133

The solution is to install the [Windows Management Framework 4.0](https://www.microsoft.com/en-us/download/details.aspx?id=40855).

## EzXML.jl cannot be installed on RedHat 6

The `zlib` library of RedHat 6, is slightly older than the library which `EzXML.jl` and `libxml2` requires.

To verify this issue, you can type in Julia

```
Libdl.dlopen(joinpath(Pkg.dir("EzXML"),"deps/usr/lib/libxml2.so"))
```

It should not return an error message. On Redhat 6.6, the following error message is returned:

```
ERROR: could not load library "/home/username/.julia/v0.6/EzXML/deps/usr/lib/libxml2.so"

/lib64/libz.so.1: version `ZLIB_1.2.3.3' not found (required by /home/divahs1/.julia/v0.6/EzXML/deps/usr/lib/libxml2.so)

Stacktrace:

 [1] dlopen(::String, ::UInt32) at ./libdl.jl:97 (repeats 2 times)
```

However, the following command should work:

```julia
 LD_LIBRARY_PATH="$HOME/.julia/v0.6/EzXML/deps/usr/lib/:$LD_LIBRARY_PATH" julia --eval  'print(Libdl.dlopen(joinpath(Pkg.dir("EzXML"),"deps/usr/lib/libxml2.so"))'
```

Lukily, EzZML.jl includes a newer version of the `zlib` library, but it does not load the library automatically.
(see also https://github.com/JuliaLang/julia/issues/7004 and https://github.com/JuliaIO/HDF5.jl/issues/97)

To make Julia use this library, a user on RedHat 6 should always start Julia with:

```bash
LD_LIBRARY_PATH="$HOME/.julia/v0.6/EzXML/deps/usr/lib/:$LD_LIBRARY_PATH" julia
```

One can also create script with the following content:

```bash
#!/bin/bash
export LD_LIBRARY_PATH="$HOME/.julia/v0.6/EzXML/deps/usr/lib/:$LD_LIBRARY_PATH"
exec /path/to/bin/julia "$@"
```

by replacing `/path/to/bin/julia` to the full path of your installation directory.
The script should be marked executable and it can be included in your Linux search [`PATH` environement variable](http://www.linfo.org/path_env_var.html). Julia can then be started by calling directly this script.

## The DIVAnd test suite fails with `automatic download failed`

Running `Pkg.test("divand")` fails with the error:

```julia
automatic download failed (error: 2147500036)
```

The test suite will download some sample data. You need to have internet access and run the test function from a directory with write access.

You can change the directory to your home directory with the julia command `cd(homedir())`.

You can check the current working directory with:

```julia
pwd()
```

## METADATA cannot be updated

`Pkg.update` fails with the error message `METADATA cannot be updated`.

If you have git installed, you can issue the command:

```bash
cd ~/.julia/v0.6/METADATA
git reset --hard
```

and then in Julia run `Pkg.update()` again.

If this does not work, then, you can also delete `~/.julia` (https://github.com/JuliaLang/julia/issues/18651#issuecomment-347579521) and in Julia enter `Pkg.init(); Pkg.update()`.


## Convert error in `divand_obs`

The full error message:

```
MethodError: Cannot `convert` an object of type divand.divand_constrain{Float32,Diagonal{Float64},SparseMatrixCSC{Float64,Int64}} to an object of type divand.divand_constrain{Float64,TR,TH} where TH<:(AbstractArray{#s370,2} where #s370<:Number) where TR<:(AbstractArray{#s371,2} where #s371<:Number)
This may have arisen from a call to the constructor divand.divand_constrain{Float64,TR,TH} where TH<:(AbstractArray{#s370,2} where #s370<:Number) where TR<:(AbstractArray{#s371,2} where #s371<:Number)(...),
since type constructors fall back to convert methods.
```

The solution is to use the same type of all input parameters: all Float32 or all Float64.



## Monthlist issue

Using comments inside list can lead to unexpected results.

This 

```julia
 monthlist = [
       [1,2,3]
       #[4,5,6]
       ]
```

should be written as

```julia
 monthlist = [
       [1,2,3]
       ]
```
