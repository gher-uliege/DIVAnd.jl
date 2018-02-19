var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "divand.jl documentation",
    "title": "divand.jl documentation",
    "category": "page",
    "text": "divand"
},

{
    "location": "index.html#divand.divandrun",
    "page": "divand.jl documentation",
    "title": "divand.divandrun",
    "category": "Function",
    "text": "divandrun(mask,pmn,xi,x,f,len,epsilon2; <keyword arguments>)\n\nPerform an n-dimensional variational analysis of the observations f located at the coordinates x. The array fi represent the interpolated field at the grid defined by the coordinates xi and the scales factors pmn.\n\nInput:\n\nmask: binary mask delimiting the domain. true is inside and false outside. For oceanographic application, this is the land-sea mask where sea is true and land is false.\npmn: scale factor of the grid. pmn is a tuple with n elements. Every      element represents the scale factor of the corresponding dimension. Its      inverse is the local resolution of the grid in a particular dimension.\nxi: tuple with n elements. Every element represents a coordinate of the final grid on which the observations are interpolated\nx: tuple with n elements. Every element represents a coordinate of the observations\nf: value of the observations minus the background estimate (m-by-1 array).   (see note)\nlen: correlation length\nepsilon2: error variance of the observations (normalized by the error variance of the background field). epsilon2 can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a difference error variance and their errors are decorrelated) or a matrix (all observations can have a difference error variance and their errors can be correlated). If epsilon2 is a scalar, it is thus the inverse of the signal-to-noise ratio.\n\nOptional input arguments specified as keyword arguments\n\nvelocity: velocity of advection constraint. The default is      no-advection constraint\nalpha: alpha is vector of coefficients multiplying various terms in the      cost function. The first element multiplies the norm.      The other i-th element of alpha multiplies the (i+1)-th derivative.      Per default, the highest derivative is m = ceil(1+neff/2) where neff is the      effective dimension of the problem (the number of dimensions with a nonzero      correlation length).\n The values of alpha is the (m+1)th row of the Pascal triangle:\n    m=0         1\n    m=1       1   1\n    m=1     1   2   1     (n=1,2)\n    m=2   1   3   3   1   (n=3,4)\n    ...\nconstraints: a structure with user specified constrain\nmoddim: modulo for cyclic dimension (vector with n elements).    Zero is used for non-cyclic dimensions. Halo points should    not be included for cyclic dimensions. For example if the first dimension    is cyclic, then the grid point corresponding to mask(1,j) should be    between mask(end,1) (left neighbor) and mask(2,j) (right neighbor)\nfracindex: fractional indices (n-by-m array). If this array is specified,    then x and xi are not used.\ninversion: direct solver (:chol for Cholesky factorization) or a    interative solver (:pcg for preconditioned conjugate gradient) can be    used.\ncompPC: function that returns a preconditioner for the primal formulation    if inversion is set to \'pcg\'. The function has the following arguments:\n     fun = compPC(iB,H,R)\nwhere iB is the inverse background error covariance, H the observation   operator and R the error covariance of the observation. The function compPC returns the   preconditioner fun(x,fx) computing fx = M  x (the inverse of M times x)   where M is a positive defined symmetric matrix [1].   Effectively, the system E⁻¹ A (E⁻¹)ᵀ (E x) = E⁻¹ b is solved for (E x) where E Eᵀ = M.   Ideally, M should this be similar to A, so that E⁻¹ A (E⁻¹)ᵀ is close to the identity matrix.\nfi0: starting field for iterative primal algorithm (same size as mask)\nf0: starting field for iterative dual algorithm (same size as the observations f)\noperatortype: Val{:sparse} for using sparse matrices (default) or Val{:MatFun} or using functions   to define the constrains.\nscale_len: true (default) if the correlation length-scale should be scaled such that the analysical   kernel reaches 0.6019072301972346 (besselk(1.,1.)) at the same distance. The kernel behaves thus similar to   the default kernel in two dimensions (alpha = [1,2,1]).\nalphabc : numerical value defining how the last grid points are stretched outward. 1, the default value mimics an infinite domain.\n\nTo have previous behaviour of finite domain use alphabc=0\n\nbtrunc : if provided defines where to truncate the calculation of the covariance matrix B. Only values up and including alpha[btrunc] will be calculated. IF the\n\n			iterative solution is calculated, the missing terms will be calculated on the fly during the conjugate gradient calulcations. Default value is none and full covariance calculation.\n\nOutput:\n\nfi: the analysed field\ns: structure with an array s.P representing the analysed error covariance\n\nNote:\n\nIf zero is not a valid first guess for your variable (as it is the case for   e.g. ocean temperature), you have to subtract the first guess from the   observations before calling divand and then add the first guess back in.\n\nExample:\n\nsee divand_simple_example.jl\n\nReferences\n\n[1]  https://en.wikipedia.org/w/index.php?title=Conjugate_gradient_method&oldid=761287292#The_preconditioned_conjugate_gradient_method\n\n\n\n"
},

{
    "location": "index.html#divand.divandgo",
    "page": "divand.jl documentation",
    "title": "divand.divandgo",
    "category": "Function",
    "text": "Compute a variational analysis of arbitrarily located observations.\n\nfi,s = divandgo(mask,pmn,xi,x,f,len,epsilon2,errormethod; ...);\n\nPerform an n-dimensional variational analysis of the observations f located at the coordinates x. The array fi represent the interpolated field at the grid defined by the coordinates xi and the scales factors pmn.\n\nInput:\n\nAs for divandrun but as a higher level routine which will automatically create windowing etc it also include the definition of the errormethod\n\nerrormethod : :cpme (clever poormans method), :none or :exact\n\nOutput:\n\nfi: the analysed field\ns: structure with an array s.P representing the analysed error covariance\n\n\n\n"
},

{
    "location": "index.html#divand.divand_averaged_bg",
    "page": "divand.jl documentation",
    "title": "divand.divand_averaged_bg",
    "category": "Function",
    "text": "fma,faanom = divand_averaged_bg(mask,pmn,xi,x,f,len,epsilon2,toaverage;moddim=[])\n\nInput:\n\nAs for divandrun, including all dimensions before averaging\n\nadditional argument:\n\ntoaverage: Array of ndims of boolean telling if in the corresponding direction averaging must be done\n\nPresently NO optional arguments from divandrun supported except moddim\n\nOutput:\n\nfma: Analysis where in the directions where toaverage is true, the same value is found\nfaanom: Data anomalies when the analysis is subtracted from the input field.\n\n\n\n"
},

{
    "location": "index.html#divand.load_mask",
    "page": "divand.jl documentation",
    "title": "divand.load_mask",
    "category": "Function",
    "text": "deprecated\n\n\n\ndeprecated\n\n\n\nxi,yi,mask = load_mask(bath_name,isglobal,xi,yi,level::Number)\n\nGenerate a land-sea mask based on the topography from the NetCDF file bathname. The parameter isglobal is true if the NetCDF file covers the whole globe and  thus the last longitude point can be considered to be right next to the first longitude point. xi and yi is a vector of the longitude and latitude grid onto which the bathymetry should be  interpolated. In the water, level is postive and in the air level is negative.\n\n\n\n"
},

{
    "location": "index.html#divand.domain",
    "page": "divand.jl documentation",
    "title": "divand.domain",
    "category": "Function",
    "text": "mask,(pm,pn),(xi,yi) = domain(bathname,bathisglobal,lonr,latr)\n\nGenerate a 2D geospatial domain based on the topography from the NetCDF file bathname.\n\n\n\nmask,(pm,pn,po),(xi,yi,zi) = domain(bathname,bathisglobal,lonr,latr,depthr)\n\nGenerate a 3D geospatial domain based on the topography from the NetCDF file bathname. if zlevel = :surface, then depthr is zero for the sea surface and positive in water (positive is down) if zlevel = :floor, then depthr is zero for the sea floor and positive in water (positive is up)\n\n\n\nmask,(pm,pn,po,pp),(xi,yi,zi,ti) = domain(bathname,bathisglobal,lonr,latr,depthr,timer)\n\nGenerate a geospatial domain based on the topography from the NetCDF file bathname.\n\n\n\n"
},

{
    "location": "index.html#divand.SDNMetadata",
    "page": "divand.jl documentation",
    "title": "divand.SDNMetadata",
    "category": "Function",
    "text": "ncglobalattrib,ncvarattrib = SDNMetadata(metadata,fi)\n\nBased on the information in the dictionary metadata and the analysed 4D field fi produce a list of NetCDF global and variable attributes for divand_save2.\n\n\n\n"
},

{
    "location": "index.html#divand.save",
    "page": "divand.jl documentation",
    "title": "divand.save",
    "category": "Function",
    "text": "save(filename,xyi,fi,varname;\n                      ncvarattrib = Dict(), ncglobalattrib = Dict(), ...)\n\nSave the result of the analysis in a NetCDF file .\n\nInput arguments\n\nfilename: the name of the NetCDF file\nxyi: tuple with n vectors. Every element in this tuple represents a coordinate of the final grid on which the observations are interpolated\nfi: the analysed field\nvarname: the name of the NetCDF variable\n\nOptional arguments:\n\nncglobalattrib: a dictionary with the global attributes\nncvarattrib: a dictionary with the variable attributes\nrelerr: relative error\n\n\n\n"
},

{
    "location": "index.html#divand.saveobs",
    "page": "divand.jl documentation",
    "title": "divand.saveobs",
    "category": "Function",
    "text": "divand.saveobs(filename,xy,ids;\n               type_save = Float32,\n               timeorigin = DateTime(1900,1,1,0,0,0),\n               )\n\nSave the location and time of the observation in the NetCDF file filename and their identifier ids. xy is a tuple with the vectors longitude, latitude, depth and time (as a vector of DateTime).\n\nOptional arguments:\n\ntype_save: the type to save the data (default Float32). However, the time  is always saved as Float64.\ntimeorigin: time origin for the time units attribute (default is\n\n1900-01-01 00:00:00)\n\n\n\n"
},

{
    "location": "index.html#divand.checkobs",
    "page": "divand.jl documentation",
    "title": "divand.checkobs",
    "category": "Function",
    "text": " checkobs(x,v,ids)\n checkobs(io::IO,x,v,ids)\n\nPrint some basic information about the coordinates x (tuple of vector) and  values v (vector) having the identifier ids (vector of strings) to check  erroneous data. It prints wheter NaNs or Infs are found and the minimum and  maximum value.\n\nIf the argument io is provided, the information is input/output stream io.\n\n\n\n"
},

{
    "location": "index.html#divand.fit_isotropic",
    "page": "divand.jl documentation",
    "title": "divand.fit_isotropic",
    "category": "Function",
    "text": "var0,len,distx,covar,fitcovar = fit_isotropic(x,v,distbin,mincount;\n                           alpha = divand.alpha_default(length(x)),\n                           len = 1.,\n                           var0 = 1.,\n                           minlen = 0.,\n                           maxlen = 10.,\n                           minvar0 = 0.,\n                           maxvar0 = 10.,\n                           tolrel = 1e-4,\n                           maxpoints = 10000,\n                           distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj))),\n                           progress = (var,len,fitness) -> nothing\n                       )\n\nDetermines the optimal correlation length len and variance (for a separation distance approaching zero) var0 of a cloud of data points with value v and coordiantes x (tuple of vectors with the coordinates).\n\nThe function can find the solution corresponding to  a local minimum which is not necessarily the global minimum.\n\nSee also empiriccovar for future information about the output parameters.\n\nOptional input parameters:\n\nalpha: if one correlation length is forced to zero during the anaylsis the values of alpha sould be set using the effective dimension. For example, if a 2D-analysis is simulated by forcing the vertical correlation length to zero, then alpha should be set to [1,2,1], otherwise alpha will be [1,3,3,1] (for for any proper 3D analysis).\nlen: initial value for the correlation length\nvar0: initial value of the variance\nminlen, maxlen: minimum and maximum value for the correlation length\nminvar0, maxvar0: minimum and maximum value for the variance\ntolrel: relative tolerance for the optimizer\nmaxpoints: maximum number of data points considered\ndistfun: function to compute the distance between point xi (vector) and   xj. Per default distun is the Eucedian distance  (xi,xj) -> sqrt(sum(abs2,xi-xj))).\nprogress: call-back function to show the progress of the optimization with  the input parameters var, len and fitness (all scalars).\n\nThe length-scale parameters and the variance have the corresponding units from  the x and v. It is therefore often necessary to provide reasonable values  for these default parameters.\n\nIf the lower bound minlen is too small, then you might get the following error:\n\nAmosException with id 4: input argument magnitude too large, complete loss of accuracy by argument reduction.\n\nIn these case, increase minlen.\n\n\n\n"
},

{
    "location": "index.html#divand.fit",
    "page": "divand.jl documentation",
    "title": "divand.fit",
    "category": "Function",
    "text": "var0opt,lensopt,distx,covar,fitcovar = fit(x,v,distbin,mincount;\n         alpha = divand.alpha_default(length(x)),\n         minlen = zeros(length(x)),\n         maxlen = ones(length(x)),\n         tolrel = 1e-4,\n         lens0 = ones(length(x)),\n         var0 = 1.,\n         minvar0 = 0.,\n         maxvar0 = 2.,\n         maxpoints = 10000,\n         distfun = (xi,xj,lens) -> sqrt(sum(abs2,(xi-xj)./lens)),\n         progress = (iter,var,len,fitness) -> nothing\n         )\n\nThe same as the function fit_isotropic except that now the correlation  length-scale lens0, minlen, maxlen, lensopt are a vectors  (one value per dimension). The distance function distfun uses an additional  parameter to compute the normalized distance.\n\nThe note of the optional parameters in divafit which also applies here.\n\n\n\n"
},

{
    "location": "index.html#divand.Anam.loglin",
    "page": "divand.jl documentation",
    "title": "divand.Anam.loglin",
    "category": "Function",
    "text": "trans,invtrans = loglin(t; epsilon = 0.)\n\nProvide the following transform log(x + epsilon) (for x < t) and its inverse. Beyond the threshold t (x ≥ t), the function is extended linearly in a  continous way.\n\ntrans,invtrans are scalar functions such that for any x (x > epsilon), x == invtrans(trans(x)).\n\nFor any array X, we have: X == invtrans.(trans.(X)).\n\n\n\n"
},

{
    "location": "index.html#divand.Anam.logit",
    "page": "divand.jl documentation",
    "title": "divand.Anam.logit",
    "category": "Function",
    "text": "trans,invtrans = logit(; min = 0., max = 1.)\n\nProvide the logit transform and its inverse. Per default the logit transform  maps values within the interval from 0 and 1. This can be changed with the  min and max parameters. Note that trans(min) = -∞ and trans(max) = +∞.  The use safety-margin might be necessary.\n\n\n\n"
},

{
    "location": "index.html#divand.divadoxml",
    "page": "divand.jl documentation",
    "title": "divand.divadoxml",
    "category": "Function",
    "text": "divand.divadoxml(filepath,varname,project,cdilist,xmlfilename;\n                 ignore_errors = false)\n\nGenerate the XML metadata file xmlfilename from the NetCDF file filepath with the  NetCDF variable varname. Project is either \"SeaDataNet\", \"EMODNET-chemistry\" or \"SeaDataCloud\". cdilist is the file from http://emodnet-chemistry.maris2.nl/download/export.zip.\n\nThe XML file contain a list of the data the originators. divadoxml will abort with an error if some combinations of EDMO code, local CDI ID are not present in the cdilist. Such errors can be ignore if ignore_errors is set to true.\n\n\n\n"
},

{
    "location": "index.html#divand.jl-documentation-1",
    "page": "divand.jl documentation",
    "title": "divand.jl documentation",
    "category": "section",
    "text": "divand.divandrun\ndivand.divandgo\ndivand.divand_averaged_bg\ndivand.load_mask\ndivand.domain\ndivand.SDNMetadata\ndivand.save\ndivand.saveobs\ndivand.checkobs\ndivand.fit_isotropic\ndivand.fit\ndivand.Anam.loglin\ndivand.Anam.logit\ndivand.divadoxml"
},

{
    "location": "index.html#Examples-1",
    "page": "divand.jl documentation",
    "title": "Examples",
    "category": "section",
    "text": "To run the example, you need to install PyPlot. In the folder examples of divand, you can run e.g. the example divand_simple_example_1D.jl by issuing:# cd(\"/path/to/divand/examples\")\ninclude(\"divand_simple_example_1D.jl\")Replace /path/to/divand/ by the installation directory of divand which is the output of Pkg.dir(\"divand\") if you installed divand using Julias package manager."
},

{
    "location": "index.html#divand.Vocab.CFVocab",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.CFVocab",
    "category": "Type",
    "text": "collection = Vocab.CFVocab(url)\n\nCreate a Dict-like object represeting the NetCDF CF Standard Name vocabulary. If the url is not provided then current CF Standard Name list http://cfconventions.org/Data/cf-standard-names/current/src/cf-standard-name-table.xml is used. Individual standard names are retirved by indexing which return an object of the type CFEntry:\n\ncollection = Vocab.CFVocab()\nentry = collection[\"sea_water_temperature\"]\n\n\n\n"
},

{
    "location": "index.html#Base.haskey-Tuple{divand.Vocab.CFVocab,Any}",
    "page": "divand.jl documentation",
    "title": "Base.haskey",
    "category": "Method",
    "text": "bool = haskey(collection::CFVocab,stdname)\n\nReturn true if stdname is part of the NetCDF CF Standard Name vocabulary collection.\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.SDNCollection",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.SDNCollection",
    "category": "Function",
    "text": "collection = SDNCollection(name)\n\nOpen the SeaDataNet collection with the name name at the URL http://www.seadatanet.org/urnurl/collection/ The collection can be indexed with brackets using the identifier.\n\nusing divand\ncollection = Vocab.SDNCollection(\"P01\")\nconcept = collection[\"PSALPR01\"]\n@show Vocab.prefLabel(concept)\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.prefLabel",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.prefLabel",
    "category": "Function",
    "text": "s = Vocab.prefLabel(c::Vocab.Concept)\n\nReturn the preferred label of a concept c\n\n\n\ns = Vocab.prefLabel(urn::AbstractString)\n\nReturn the preferred label of a concept usings it URN (Uniform Resource Name)\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.altLabel",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.altLabel",
    "category": "Function",
    "text": "s = Vocab.altLabel(c::Vocab.Concept)\n\nReturn the alternative label of a concept c\n\n\n\ns = Vocab.altLabel(urn::AbstractString)\n\nReturn the alternative label of a concept usings it URN (Uniform Resource Name)\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.notation",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.notation",
    "category": "Function",
    "text": "s = Vocab.notation(c::Vocab.Concept)\n\nReturn the identifier of a concept c\n\n\n\ns = Vocab.notation(urn::AbstractString)\n\nReturn the identifier of a concept usings it URN (Uniform Resource Name)\n\n\n\n"
},

{
    "location": "index.html#Base.find-Tuple{divand.Vocab.Concept,Any,Any}",
    "page": "divand.jl documentation",
    "title": "Base.find",
    "category": "Method",
    "text": "find(c::Concept,name,collection)\n\nReturn a list of related concepts in the collection collection. name can be the string \"related\", \"narrower\", \"broader\".\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.description",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.description",
    "category": "Function",
    "text": "str = description(entry::CFEntry)\nstr = canonical_units(entry::CFEntry)\n\nReturn the description or the canonical units of the `entry`.\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.canonical_units",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.canonical_units",
    "category": "Function",
    "text": "str = description(entry::CFEntry)\nstr = canonical_units(entry::CFEntry)\n\nReturn the description or the canonical units of the `entry`.\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.splitURL",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.splitURL",
    "category": "Function",
    "text": "collection,tag,key = Vocab.splitURL(url)\n\nSplit a concept URL into collection, tag and key. url must finishe with a slash.\n\n\n\n"
},

{
    "location": "index.html#Vocabulary-1",
    "page": "divand.jl documentation",
    "title": "Vocabulary",
    "category": "section",
    "text": "urn_strVocab.CFVocab\nhaskey(collection::Vocab.CFVocab,stdname)\nVocab.SDNCollection\nVocab.prefLabel\nVocab.altLabel\nVocab.notation\nVocab.find(c::Vocab.Concept,name,collection)\nVocab.description\nVocab.canonical_units\nVocab.splitURL"
},

{
    "location": "index.html#Information-for-developpers-1",
    "page": "divand.jl documentation",
    "title": "Information for developpers",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Update-the-documentation-1",
    "page": "divand.jl documentation",
    "title": "Update the documentation",
    "category": "section",
    "text": "InstallPkg.add(\"Documenter\")"
},

]}
