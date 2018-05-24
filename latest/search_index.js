var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "divand.jl documentation",
    "title": "divand.jl documentation",
    "category": "page",
    "text": "divand"
},

{
    "location": "index.html#divand.diva3d",
    "page": "divand.jl documentation",
    "title": "divand.diva3d",
    "category": "function",
    "text": "residuals = diva3d(xi,x,value,len,epsilon2,filename,varname)\n\nCreate a 3D analysis (or a series of 3D analyses) with DIVAnd using the observations value (vector) at the locations x (tuple of vectors) onto the regular grid defined by the vectors xi using the scaled observational error variance  epsilon2 and the correlation length len. The result will be saved in the NetCDF file filename under the variable varname.\n\nInputs\n\nxi: tuple with n elements. Every element represents a coordinate of the final grid on which the observations are interpolated\nx: tuple with n elements. Every element represents a coordinate of the observations\nvalue: value of the observations\nlen: tuple with n elements. Every element represents the correlation length.  If fitcorrlen is true, then len can be the empty tuple () or a tuple containing  3 arrays of normalized correlation lengths which will be multiplied by the  horizontal and vertical correlation lengths.\nepsilon2: error variance of the observations (normalized by the error variance of the background field). epsilon2 can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a different error variance and their errors are decorrelated) or a matrix (all observations can have a different error variance and their errors can be correlated). If epsilon2 is a scalar, it is thus the inverse of the signal-to-noise ratio.\nfilename: The output NetCDF filename.\nvarname: The name of the variable (used in the NetCDF file).\n\nOptional input arguments:\n\nbathname: path to the NetCDF bathymetry (default ../../divand-example-data/Global/Bathymetry/gebco_30sec_16.nc relative to this source file)\nbathisglobal: true (default) is the bathymetry is a global data set\nplotres: Call-back routine for plotting ((timeindex,sel,fit,erri) -> nothing)\ntimeorigin: Time origin (default DateTime(1900,1,1,0,0,0))\nmoddim: modulo for cyclic dimension (vector with n elements).    Zero is used for non-cyclic dimensions. Halo points should    not be included for cyclic dimensions. For example if the first dimension    is cyclic, then the grid point corresponding to mask[1,j] should be    between mask[end,1] (left neighbor) and mask[2,j] (right neighbor). The default is [0,0,0],\nzlevel: :surface (default) for surface analysis and :floor for analysis from the bottom floor.\nncvarattrib: dictionary of NetCDF variable attributes.\nncglobalattrib: dictionary of NetCDF global attributes.\ntransform: Anamorphosis transformation function (default: Anam.notransform()).\nfitcorrlen: true of the correlation length is determined from the observation (default false).\nfithorz_param: dictionary with additional optional parameters for fithorzlen.\nfitvert_param: dictionary with additional optional parameters for fitvertlen.\ndistfun: function to compute the distance (default (xi,xj) -> divand.distance(xi[2],xi[1],xj[2],xj[1])).\nmask: if different from nothing, then this mask overrides land-sea mask based on the bathymetry\n\n(default nothing).\n\nbackground: if different from nothing, then this parameter allows one\n\nto load the background from a call-back function (default nothing).\n\nbackground_espilon2_factor: multiplication for epsilon2 when computing the background (default 10.).\nmemtofit: keyword controlling how to cut the domain depending on the memory   remaining available for inversion. It is not total memory (default 3).\nniter_e: Number of iterations to estimate the optimal scale factor of  epsilon2 using Desroziers et al. 2005 (doi: 10.1256/qj.05.108). The default   is 1 (i.e. no optimization is done).\n\nAny additional keywoard arguments understood by divandgo can also be used here (e.g. velocity constrain)\n\n\n\n"
},

{
    "location": "index.html#divand.divandrun",
    "page": "divand.jl documentation",
    "title": "divand.divandrun",
    "category": "function",
    "text": "divandrun(mask,pmn,xi,x,f,len,epsilon2; <keyword arguments>)\n\nPerform an n-dimensional variational analysis of the observations f located at the coordinates x. The array fi represent the interpolated field at the grid defined by the coordinates xi and the scales factors pmn.\n\nInput:\n\nmask: binary mask delimiting the domain. true is inside and false outside.\n\nFor oceanographic application, this is the land-sea mask where sea is true and land is false.\n\npmn: scale factor of the grid. pmn is a tuple with n elements. Every  element represents the scale factor of the corresponding dimension. Its  inverse is the local resolution of the grid in a particular dimension.  For example, in two dimensions, pmn is a tuple (pm,pn) where pm is  the inverse of the local resolution in first dimension and pn is the the inverse  of the local resolution in second dimension.\nxi: tuple with n elements. Every element represents a coordinate of the final grid on which the observations are interpolated.\nx: tuple with n elements. Every element represents a coordinate of the observations.\nf: value of the observations minus the background estimate (vector of m elements where m is the number of observations). See also note.\nlen: tuple with n elements. Every element represents the correlation length for a given dimension.\nepsilon2: error variance of the observations (normalized by the error variance of the background field). epsilon2 can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a different error variance and their errors are decorrelated) or a matrix (all observations can have a different error variance and their errors can be correlated). If epsilon2 is a scalar, it is thus the inverse of the signal-to-noise ratio.\n\nOptional input arguments specified as keyword arguments\n\nvelocity: velocity of advection constraint. The default is      no-advection constraint\nalpha: alpha is vector of coefficients multiplying various terms in the      cost function. The first element multiplies the norm.      The other i-th element of alpha multiplies the (i+1)-th derivative.      Per default, the highest derivative is m = ceil(1+neff/2) where neff is the      effective dimension of the problem (the number of dimensions with a nonzero      correlation length) and ceil is the ceiling function (rounding up).\n\n   The values of alpha is the (m+1)th row of the Pascal triangle:\n      m=0         1\n      m=1       1   1\n      m=1     1   2   1     (n=1,2)\n      m=2   1   3   3   1   (n=3,4)\n      ...\n\nconstraints: a structure with user specified constraints (see divand_addc).\nmoddim: modulo for cyclic dimension (vector with n elements).    Zero is used for non-cyclic dimensions. One should not include a boundary    zone (sometimes called a ghost zone or halo) for cyclic dimensions.    For example if the first dimension    is cyclic, then the grid point corresponding to mask[1,j] should be    between mask[end,1] (left neighbor) and mask[2,j] (right neighbor).\nfracindex: fractional indices (n-by-m array). If this array is specified,    then x and xi are not used.\ninversion: direct solver (:chol for Cholesky factorization) or an    interative solver (:pcg for preconditioned conjugate gradient [1]) can be    used.\ncompPC: function that returns a preconditioner for the primal formulation    if inversion is set to \'pcg\'. The function has the following arguments:\n     fun = compPC(iB,H,R)\nwhere iB is the inverse background error covariance, H the observation   operator and R the error covariance of the observation. The function compPC returns the   preconditioner fun(x,fx) computing fx = M  x (the inverse of M times x)   where M is a positive defined symmetric matrix [1].   Effectively, the system E⁻¹ A (E⁻¹)ᵀ (E x) = E⁻¹ b is solved for (E x) where E Eᵀ = M.   Ideally, M should this be similar to A, so that E⁻¹ A (E⁻¹)ᵀ is close to the identity matrix.\nfi0: starting field for iterative primal algorithm (same size as mask).\nf0: starting field for iterative dual algorithm (same size as the observations f).\noperatortype: Val{:sparse} for using sparse matrices (default) or Val{:MatFun} or using functions   to define the constrains.\nscale_len: true (default) if the correlation length-scale should be scaled   such that the analytical   kernel reaches 0.6019072301972346 (besselk(1.,1.)) at the same distance   than in 2D. The kernel behaves thus similar to   the default kernel in two dimensions (alpha = [1,2,1]).\nalphabc : numerical value defining how the last grid points are stretched outward.  If alphabc is 1, the default value mimics an infinite domain.  To have previous behaviour of finite domain use alphabc equal to 0.\nbtrunc : if provided defines where to truncate the calculation of the   covariance matrix B. Only values up and including alpha[btrunc] will be   calculated. If the iterative solution is calculated, the missing terms will   be calculated on the fly during the conjugate gradient calulcations. Default value is none and full covariance calculation.\n\nOutput:\n\nfi: the analysed field\ns: structure with an array s.P representing the analysed error covariance\n\nNote:\n\nIf zero is not a valid first guess for your variable (as it is the case for   e.g. ocean temperature), you have to subtract the first guess from the   observations before calling divand and then add the first guess back in.\n\nExample:\n\nsee divand_simple_example.jl\n\nReferences\n\n[1]  https://en.wikipedia.org/w/index.php?title=Conjugate_gradient_method&oldid=761287292#The_preconditioned_conjugate_gradient_method\n\n\n\n"
},

{
    "location": "index.html#divand.divandgo",
    "page": "divand.jl documentation",
    "title": "divand.divandgo",
    "category": "function",
    "text": "fi, erri, residuals, qcvalues, scalefactore = divandgo(mask,pmn,xi,x,f,len,epsilon2,errormethod; ...);\n\nInput:\n\nSame arguments as divandrun with in addition\nerrormethod :   you have the choice between :cpme (clever poorman\'s method, default method if parameter not provided), :none or :exact (only available if windowed analysis are done with divandrun)\nMEMTOFIT=: keyword controlling how to cut the domain depending on the memory remaining available for inversion (not total memory)\nRTIMESONESCALES= : if you provide a tuple of length scales, data are weighted differently depending on the numbers of neighbours they have. See weight_RtimesOne for details \nQCMETHOD= : if you provide a qc method parameter, quality flags are calculated. See divand_cv for details\n\nOutput:\n\nfi: the analysed field\nerri: relative error field on the same grid as fi. () if errormethod is fixed to :none\nresiduals: array of residuals at data points. For points not on the grid or on land: NaN\nqcvalues: if QCMETHOD= is provided, the output array contains the quality flags otherwise qcvalues is (). For points on land or not on the grid: 0\nscalefactore: Desroziers et al. 2005 (doi: 10.1256/qj.05.108) scale factor for epsilon2\n\nPerform an n-dimensional variational analysis of the observations f located at the coordinates x. The array fi represent the interpolated field at the grid defined by the coordinates xi and the scales factors pmn. \n\nIMPORTANT: divandgo is very similar to divandrun and is only interesting to use if divandrun cannot fit into memory or if you want to parallelize. (In the latter case do not forget to define the number of workers; see addprocs for example)\n\n\n\n"
},

{
    "location": "index.html#divand.divand_averaged_bg",
    "page": "divand.jl documentation",
    "title": "divand.divand_averaged_bg",
    "category": "function",
    "text": "fma,faanom = divand_averaged_bg(mask,pmn,xi,x,f,len,epsilon2,toaverage;moddim=[])\n\nInput:\n\nAs for divandrun, including all dimensions before averaging\n\nadditional argument:\n\ntoaverage: Array of ndims of boolean telling if in the corresponding direction averaging must be done\n\nPresently NO optional arguments from divandrun supported except moddim\n\nOutput:\n\nfma: Analysis where in the directions where toaverage is true, the same value is found\nfaanom: Data anomalies when the analysis is subtracted from the input field.\n\n\n\n"
},

{
    "location": "index.html#divand.SDNMetadata",
    "page": "divand.jl documentation",
    "title": "divand.SDNMetadata",
    "category": "function",
    "text": "ncglobalattrib,ncvarattrib = SDNMetadata(metadata,fi)\n\nBased on the information in the dictionary metadata and the analysed 4D field fi produce a list of NetCDF global and variable attributes for divand_save2.\n\n\n\n"
},

{
    "location": "index.html#divand.save",
    "page": "divand.jl documentation",
    "title": "divand.save",
    "category": "function",
    "text": "save(filename,xyi,fi,varname;\n                      ncvarattrib = Dict(), ncglobalattrib = Dict(), ...)\n\nSave the result of the analysis in a NetCDF file .\n\nInput arguments\n\nfilename: the name of the NetCDF file\nxyi: tuple with n vectors. Every element in this tuple represents a coordinate of the final grid on which the observations are interpolated\nfi: the analysed field\nvarname: the name of the NetCDF variable\n\nOptional arguments:\n\nncglobalattrib: a dictionary with the global attributes\nncvarattrib: a dictionary with the variable attributes\nrelerr: relative error\n\n\n\n"
},

{
    "location": "index.html#divand.loadbigfile",
    "page": "divand.jl documentation",
    "title": "divand.loadbigfile",
    "category": "function",
    "text": "value,lon,lat,depth,time,obsid = loadbigfile(filename)\n\nLoad data from the text file filename and returns vectors with the value, longitude, latitude, depth and time (as DateTime). A list string identifiers is also returned.\n\n\n\n"
},

{
    "location": "index.html#divand.checkobs",
    "page": "divand.jl documentation",
    "title": "divand.checkobs",
    "category": "function",
    "text": " checkobs(x,v,ids)\n checkobs(io::IO,x,v,ids)\n\nPrint some basic information about the coordinates x (tuple of vector) and values v (vector) having the identifier ids (vector of strings) to check erroneous data. It prints wheter NaNs or Infs are found and the minimum and maximum value.\n\nIf the argument io is provided, the information is input/output stream io.\n\n\n\n"
},

{
    "location": "index.html#divand.smoothfilter",
    "page": "divand.jl documentation",
    "title": "divand.smoothfilter",
    "category": "function",
    "text": "ff = smoothfilter(x,f,scale)\n\nSmooth the function f defined on x by solving the diffusion equation\n\n∂ₜ ϕ = ν ∂²ₓ ϕ\n\nscale is the spatial scales of the removed length-scales. It is defined as 2Tν  where T is the integration time.\n\nIt uses the Greens functions for 1D diffusion: 1/sqrt(4 π ν t) * exp(-x^2 / (4νt))\n\n\n\n"
},

{
    "location": "index.html#divand.Anam.loglin",
    "page": "divand.jl documentation",
    "title": "divand.Anam.loglin",
    "category": "function",
    "text": "trans,invtrans = loglin(t; epsilon = 0.)\n\nProvide the following transform log(x + epsilon) (for x < t) and its inverse. Beyond the threshold t (x ≥ t), the function is extended linearly in a  continous way.\n\ntrans,invtrans are scalar functions such that for any x (x > epsilon), x == invtrans(trans(x)).\n\nFor any array X, we have: X == invtrans.(trans.(X)).\n\n\n\n"
},

{
    "location": "index.html#divand.Anam.logit",
    "page": "divand.jl documentation",
    "title": "divand.Anam.logit",
    "category": "function",
    "text": "trans,invtrans = logit(; min = 0., max = 1.)\n\nProvide the logit transform and its inverse. Per default the logit transform  maps values within the interval from 0 and 1. This can be changed with the  min and max parameters. Note that trans(min) = -∞ and trans(max) = +∞.  The use safety-margin might be necessary.\n\n\n\n"
},

{
    "location": "index.html#divand.divadoxml",
    "page": "divand.jl documentation",
    "title": "divand.divadoxml",
    "category": "function",
    "text": "divand.divadoxml(filepath,varname,project,cdilist,xmlfilename;\n                 ignore_errors = false)\n\nGenerate the XML metadata file xmlfilename from the NetCDF file filepath with the  NetCDF variable varname. Project is either \"SeaDataNet\", \"EMODNET-chemistry\" or \"SeaDataCloud\". cdilist is the file from http://emodnet-chemistry.maris2.nl/download/export.zip.\n\nThe XML file contain a list of the data the originators. divadoxml will abort with an error if some combinations of EDMO code, local CDI ID are not present in the cdilist. Such errors can be ignore if ignore_errors is set to true.\n\n\n\n"
},

{
    "location": "index.html#divand.random",
    "page": "divand.jl documentation",
    "title": "divand.random",
    "category": "function",
    "text": "field = divand.random(mask,pmn,len,Nens)\n\nCreate Nens random fields with the correlation length len in  a domain with the mask mask and the metric pmn.\n\nSee divand.divandrun for more information about these parameters.\n\n\n\n"
},

{
    "location": "index.html#divand.distance",
    "page": "divand.jl documentation",
    "title": "divand.distance",
    "category": "function",
    "text": "d = distance(lat1,lon1,lat2,lon2)\n\nCompute the great-circle distance between the points (lat1,lon1) and (lat2,lon2). The units of all input and output parameters are degrees.\n\n\n\nd = distance([lon1,lat1],[lon2,lat2])\n\nThe same as distance(lat1,lon1,lat2,lon2) but there the arguments are vectors  and the order is longitude then latitude.\n\nThe units of all input and output parameters are degrees.\n\n\n\n"
},

{
    "location": "index.html#divand.interp",
    "page": "divand.jl documentation",
    "title": "divand.interp",
    "category": "function",
    "text": "f = interp(xi,fi,x)\n\nInterpolate field fi (n-dimensional array) defined at xi (tuble of n-dimensional arrays or vectors) onto grid x (tuble of n-dimensional arrays). The grid in xi must be align with the axis (e.g. produced by divand.ndgrid).\n\n\n\n"
},

{
    "location": "index.html#divand.backgroundfile",
    "page": "divand.jl documentation",
    "title": "divand.backgroundfile",
    "category": "function",
    "text": "fun = backgroundfile(fname,varname)\n\nReturn a function fun which is used in DIVAnd to make  anomalies out of observations based relative to the field defined in the NetCDF variable varname in the NetCDF file  fname. It is assumed that the NetCDF variables has the variable lon, lat and depth. And that the NetCDF variable is defined on the  same grid as the analysis.\n\n\n\n"
},

{
    "location": "index.html#divand.Quadtrees.checkduplicates",
    "page": "divand.jl documentation",
    "title": "divand.Quadtrees.checkduplicates",
    "category": "function",
    "text": "dupl = checkduplicates(x,value,delta,deltavalue)\n\nBased the coordinates x (a tuple of longitude lons, latitudes lats, depth (zs)  and time (times vector of DateTime)) check of points who are in the same spatio-temporal bounding  box of a length delta. delta is a vector with 4 elements corresponding to  longitude, latitude, depth and time (in days). dupl a vector of vectors containing indices of the duplicates.\n\n\n\ndupl = checkduplicates(x1,value1,x2,v2,value2,delta,deltavalue)\n\nReport duplicate of observation in data set (x2,v2) which are also in data set  (x1,v1). x1 and x2 is a tuple of vectors with the cooridantes and v1 and v2 the  corresponding values. \n\n\n\n"
},

{
    "location": "index.html#divand.jl-documentation-1",
    "page": "divand.jl documentation",
    "title": "divand.jl documentation",
    "category": "section",
    "text": "divand.diva3d\ndivand.divandrun\ndivand.divandgo\ndivand.divand_averaged_bg\ndivand.SDNMetadata\ndivand.save\ndivand.loadbigfile\ndivand.checkobs\ndivand.smoothfilter\ndivand.Anam.loglin\ndivand.Anam.logit\ndivand.divadoxml\ndivand.random\ndivand.distance\ndivand.interp\ndivand.backgroundfile\ndivand.Quadtrees.checkduplicates"
},

{
    "location": "index.html#divand.load_bath",
    "page": "divand.jl documentation",
    "title": "divand.load_bath",
    "category": "function",
    "text": "xi,yi,bath = divand.load_bath(bath_name,isglobal,xi,yi)\n\nLoad the bathymetry from the NetCDF file bathname. The parameter isglobal is true if the NetCDF file covers the whole globe and thus the last longitude point can be considered to be right next to the first longitude point. xi and yi are vectors containing the longitude and latitude grid onto which the bathymetry should be interpolated.\n\n\n\n"
},

{
    "location": "index.html#divand.extract_bath",
    "page": "divand.jl documentation",
    "title": "divand.extract_bath",
    "category": "function",
    "text": "bx,by,b = divand.extract_bath(bath_name,isglobal,xi,yi)\n\nExtract the bathymetry from the NetCDF file bathname. The parameter isglobal  is true if the NetCDF file covers the whole globe and thus the last longitude point can be considered to be right next to the first longitude point. xi and yi are vectors defining the bounding box of the data. No interpolation is performed.\n\nConvention: b is positive in the water and negative in the air.\n\n\n\n"
},

{
    "location": "index.html#divand.load_mask",
    "page": "divand.jl documentation",
    "title": "divand.load_mask",
    "category": "function",
    "text": "xi,yi,mask = load_mask(bath_name,isglobal,xi,yi,level::Number)\n\nGenerate a land-sea mask based on the topography from the NetCDF file bathname. The parameter isglobal is true if the NetCDF file covers the whole globe and thus the last longitude point can be considered to be right next to the first longitude point. xi and yi are vectors containing the longitude and latitude grid onto which the bathymetry should be interpolated.\n\nConvention: in the water, level is positive and in the air level is negative.\n\n\n\n"
},

{
    "location": "index.html#divand.divand_metric",
    "page": "divand.jl documentation",
    "title": "divand.divand_metric",
    "category": "function",
    "text": "pm,pn = divand_metric(lon,lat)\n\nCompute metric scale factors pm and pn based on the arrays longitude lon and latitude lat. The variables pm and pn represent the inverse of the local resolution in meters using the mean Earth radius.\n\n\n\n"
},

{
    "location": "index.html#divand.domain",
    "page": "divand.jl documentation",
    "title": "divand.domain",
    "category": "function",
    "text": "mask,(pm,pn),(xi,yi) = domain(bathname,bathisglobal,lonr,latr)\n\nGenerate a 2D geospatial domain based on the topography from the NetCDF file bathname.\n\n\n\nmask,(pm,pn,po),(xi,yi,zi) = domain(bathname,bathisglobal,lonr,latr,depthr)\n\nGenerate a 3D geospatial domain based on the topography from the NetCDF file bathname. If zlevel is :surface, then depthr is zero for the sea surface and  positive in water (positive is down). If zlevel is :floor, then depthr is  zero for the sea floor and positive in water (positive is up)\n\n\n\nmask,(pm,pn,po,pp),(xi,yi,zi,ti) = domain(bathname,bathisglobal,lonr,latr,depthr,timer)\n\nGenerate a geospatial domain based on the topography from the NetCDF file bathname.\n\n\n\n"
},

{
    "location": "index.html#divand.divand_rectdom",
    "page": "divand.jl documentation",
    "title": "divand.divand_rectdom",
    "category": "function",
    "text": "mask,pmn,xyi = divand_rectdom(coord1,coord2,...)\n\nCreate a \"rectangular\" domain in n dimensions with the coordinates coord1 coord2... assuming a Catersian metric. This functions returns the mask mask, the coordinates (xi,yi,...) and the metric (pm,pn...).\n\nFor example:\n\njulia> mask,(pm,pn),(xi,yi) = divand_rectdom(linspace(0,1,50),linspace(0,1,50))\n\n\n\n"
},

{
    "location": "index.html#divand.divand_squaredom",
    "page": "divand.jl documentation",
    "title": "divand.divand_squaredom",
    "category": "function",
    "text": "mask,pmn,xyi = divand_squaredom(n,coord)\n\nCreate a \"square\" domain in n dimensions with the coordinates coord assuming a Catersian metric. This functions returns the mask mask, the coordinates (xi,yi,...) and the metric (pm,pn...).\n\nExample\n\nmask,(pm,pn),(xi,yi) = divand_squaredom(2,linspace(0,1,50))\n\n\n\n"
},

{
    "location": "index.html#divand.TimeSelectorYW",
    "page": "divand.jl documentation",
    "title": "divand.TimeSelectorYW",
    "category": "function",
    "text": "TS = TimeSelectorYW(years,yearwindow,monthlists)\n\nThe structure TS handles the time aggregation based on years and monthlists. It is similar to TimeSelectorYearListMonthList except that the elements of yearlists are centred around years and span yearwindow years. yearlists is in fact constructed by adding and subtracting yearwindow/2 to every element of years.\n\n\n\n"
},

{
    "location": "index.html#divand.TimeSelectorYearListMonthList",
    "page": "divand.jl documentation",
    "title": "divand.TimeSelectorYearListMonthList",
    "category": "type",
    "text": "TS = TimeSelectorYearListMonthList(yearlists,monthlists)\n\nThe structure TS handles the time aggregation based on yearlists and monthlists. yearlists is a vector of ranges (containing start and end years), for example [1980:1990,1990:2000,2000:2010].\n\nmonthlists is a vector of two-element vector (containing start and end months), for example [1:3,4:6,7:9,10:12]\n\nIf a month range spans beyond December, then all Months must be specified, e.g. example [2:4,5:6,7:9,[10,11,12,1]] or [2:4,5:6,7:9,[10:12;1]]. However using [2:4,5:6,7:9,10:1] (bug!) will result in an empty month range.\n\nExample\n\nseasonal climatology using all data from 1900 to 2017\n\nfor winter (December-February), spring, summer, autumn\n\nTS = divand.TimeSelectorYearListMonthList([1900:2017],[[12,1,2],[3,4,5],[6,7,8],[9,10,11]])\n\n\n\n"
},

{
    "location": "index.html#Bathymetry-and-spatial-temporal-domain-1",
    "page": "divand.jl documentation",
    "title": "Bathymetry and spatial-temporal domain",
    "category": "section",
    "text": "divand.load_bath\ndivand.extract_bath\ndivand.load_mask\ndivand.divand_metric\ndivand.domain\ndivand.divand_rectdom\ndivand.divand_squaredom\ndivand.TimeSelectorYW\ndivand.TimeSelectorYearListMonthList"
},

{
    "location": "index.html#divand.saveobs",
    "page": "divand.jl documentation",
    "title": "divand.saveobs",
    "category": "function",
    "text": "divand.saveobs(filename,xy,ids;\n               type_save = Float32,\n               timeorigin = DateTime(1900,1,1,0,0,0),\n               )\n\nSave the location and time of the observation in the NetCDF file filename and their identifier ids. xy is a tuple with the vectors longitude, latitude, depth and time (as a vector of DateTime).\n\nOptional arguments:\n\ntype_save: the type to save the data (default Float32). However, the time  is always saved as Float64.\ntimeorigin: time origin for the time units attribute (default is\n\n1900-01-01 00:00:00)\n\n\n\ndivand.saveobs(filename,varname,value,xy,ids;\n               type_save = Float32,\n               timeorigin = DateTime(1900,1,1,0,0,0),\n               )\n\nSave value and the location and time of the observation in the NetCDF file filename and their identifier ids. xy is a tuple with the vectors longitude, latitude, depth and time (as a vector of DateTime). The values will be saved in the  variable called varname.\n\nOptional arguments:\n\ntype_save: the type to save the data (default Float32). However, the time  is always saved as Float64.\ntimeorigin: time origin for the time units attribute (default is\n\n1900-01-01 00:00:00)\n\n\n\n"
},

{
    "location": "index.html#divand.loadobs",
    "page": "divand.jl documentation",
    "title": "divand.loadobs",
    "category": "function",
    "text": "value,lon,lat,depth,time,obsid = loadobs(T,filename,varname)\n\nLoad the variable varname from the NetCDF file filename. Coordinates (the NetCDF variables \"obslon\", \"obslat\", \"obsdepth\"), time (\"obstime\") and identifies (\"obsids\") will also be loaded. Numeric output arguments will have the type T.\n\n\n\n"
},

{
    "location": "index.html#divand.NCSDN.load",
    "page": "divand.jl documentation",
    "title": "divand.NCSDN.load",
    "category": "function",
    "text": "data,lon,lat,z,time,ids = load(T,fname::TS,param; qualityflags = [GOOD_VALUE, PROBABLY_GOOD_VALUE]) where TS <: AbstractString\n\n\n\ndata,lon,lat,z,time,ids = SDN.load(T,fnames,param; qualityflags = ...)\n\nLoad all data in the vector of file names fnames corresponding to the parameter  param as the data type T. Only the data with the quality flags  SDN.good_data and SDN.probably_good_data are loaded per default. The output parameters correspondata to the data, longitude, latitude, depth, time (as DateTime) and an identifier (as String).\n\n\n\n"
},

{
    "location": "index.html#divand.NCSDN.loadvar",
    "page": "divand.jl documentation",
    "title": "divand.NCSDN.loadvar",
    "category": "function",
    "text": "data = loadvar(ds,param;\n               fillvalue::T = NaN,\n               qualityflags = [GOOD_VALUE, PROBABLY_GOOD_VALUE],\n               qfname = param * QC_SUFFIX,\n               )\n\nLoad the NetCDF variable param from the NCDataset ds.  Data points not having the provide quality flags will be masked by fillvalue. qfname is the NetCDF variable name for the quality flags.\n\n\n\n"
},

{
    "location": "index.html#divand.ODVspreadsheet.loaddata",
    "page": "divand.jl documentation",
    "title": "divand.ODVspreadsheet.loaddata",
    "category": "function",
    "text": "data = loaddata(sheet,profile,locname,fillvalue; fillmode = :repeat)\n\nLoad a single column referred by the local name locname in the profile profile from the ODV spreadsheet sheet. Empty values are either replaced by fillvalue (if fillmode is :fill) or the previous value if repeated (if fillmode is :repeat)\n\n\n\n"
},

{
    "location": "index.html#divand.ODVspreadsheet.parsejd",
    "page": "divand.jl documentation",
    "title": "divand.ODVspreadsheet.parsejd",
    "category": "function",
    "text": "dt = parsejd(t)\n\nConvert a Chronological Julian Day Number to a DateTime object. The reference value is taken from Chronological Julian Date\n\nFrom the SDN standard: \"A real number representing the Chronological Julian Date, which is defined as the time elapsed in days from 00:00 on January 1 st 4713 BC. ... \"\n\nThe time origin is _not_ noon (12:00) on Monday, January 1, 4713 BC as for the Julia Date Number.\n\n\n\n"
},

{
    "location": "index.html#divand.ODVspreadsheet.myparse",
    "page": "divand.jl documentation",
    "title": "divand.ODVspreadsheet.myparse",
    "category": "function",
    "text": "v = myparse(T,s)\n\nParse the string s as a type T. Unlike Julia\'s parse function an error message contains the string s (which could not be parsed) for debugging.\n\n\n\n"
},

{
    "location": "index.html#Load-observations-1",
    "page": "divand.jl documentation",
    "title": "Load observations",
    "category": "section",
    "text": "divand.saveobs\ndivand.loadobs\ndivand.NCSDN.load\ndivand.NCSDN.loadvar\ndivand.ODVspreadsheet.loaddata\ndivand.ODVspreadsheet.parsejd\ndivand.ODVspreadsheet.myparse"
},

{
    "location": "index.html#divand.fit_isotropic",
    "page": "divand.jl documentation",
    "title": "divand.fit_isotropic",
    "category": "function",
    "text": "var0,len,distx,covar,fitcovar = fit_isotropic(x,v,distbin,mincount;\n                           alpha = divand.alpha_default(length(x)),\n                           minlen = 0.,\n                           maxlen = 10.,\n                           tolrel = 1e-4,\n                           maxpoints = 10000,\n                           nmean = 100,\n                           distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj))),\n                           progress = (iter,var,len,fitness) -> nothing\n                       )\n\nDetermines the optimal correlation length len and variance (for a separation distance approaching zero) var0 of a cloud of data points with value v and coordiantes x (tuple of vectors with the coordinates).\n\nThe function can find the solution corresponding to a local minimum which is not necessarily the global minimum.\n\nSee also empiriccovar for future information about the output parameters.\n\nOptional input parameters:\n\nalpha: if one correlation length is forced to zero during the anaylsis the values of alpha sould be set using the effective dimension. For example, if a 2D-analysis is simulated by forcing the vertical correlation length to zero, then alpha should be set to [1,2,1], otherwise alpha will be [1,3,3,1] (for any proper 3D analysis).\nlen: initial value for the correlation length.\nminlen, maxlen: minimum and maximum values for the correlation length.\ntolrel: relative tolerance for the optimizer.\nmaxpoints: maximum number of data points considered.\nnmean: the number of times an empirical covariance is estimated.  The average covariance is used for the fitting.\ndistfun: function to compute the distance between point xi (vector) and  xj. Per default distfun is the Euclidian distance: (xi,xj) -> sqrt(sum(abs2,xi-xj))).\nprogress: call-back function to show the progress of the optimization with the input parameters iter, var, len and fitness (all scalars).\n\nThe length-scale parameters and the variance have the corresponding units from the x and v. It is therefore often necessary to provide reasonable values for these default parameters.\n\nThe algorithm used to estimate the correlation-length and variance is based on randomly choosen points. Therefore the result can be different if the function is invoked repeately. If nmean is increased, then these statistical fluctuations should decrease (for a not too large value of mincount, i.e. about 100 for most cases).\n\nIf the lower bound minlen is too small, then you might get the following error:\n\nAmosException with id 4: input argument magnitude too large, complete loss of accuracy by argument reduction.\n\nIn these case, increase minlen.\n\n\n\n"
},

{
    "location": "index.html#divand.fit",
    "page": "divand.jl documentation",
    "title": "divand.fit",
    "category": "function",
    "text": "var0opt,lensopt,distx,covar,fitcovar = fit(x,v,distbin,mincount;\n         alpha = divand.alpha_default(length(x)),\n         minlen = zeros(length(x)),\n         maxlen = ones(length(x)),\n         tolrel = 1e-4,\n         lens0 = ones(length(x)),\n         var0 = 1.,\n         minvar0 = 0.,\n         maxvar0 = 2.,\n         maxpoints = 10000,\n         distfun = (xi,xj,lens) -> sqrt(sum(abs2,(xi-xj)./lens)),\n         progress = (iter,var,len,fitness) -> nothing\n         )\n\nThe same as the function fit_isotropic except that now the correlation length-scale lens0, minlen, maxlen, lensopt are vectors (one value per dimension). The distance function distfun uses an additional parameter to compute the normalized distance.\n\nThe note of the optional parameters in divafit also applies here.\n\n\n\n"
},

{
    "location": "index.html#divand.divand_cv",
    "page": "divand.jl documentation",
    "title": "divand.divand_cv",
    "category": "function",
    "text": "bestfactorl,bestfactore, cvval,cvvalues, x2Ddata,y2Ddata,cvinter,xi2D,yi2D = divand_cv(mask,pmn,xi,x,f,len,epsilon2,nl,ne,method;...);\n\nInput\n\nSame as for divandrun with three more parameters nl,ne and method\n\nmask: binary mask delimiting the domain. true is inside and false outside. For oceanographic application, this is the land-sea mask.\npmn: scale factor of the grid. pmn is a tuple with n elements. Every      element represents the scale factor of the corresponding dimension. Its      inverse is the local resolution of the grid in a particular dimension.\nxi: tuple with n elements. Every element represents a coordinate of the final grid on which the observations are interpolated\nx: tuple with n elements. Every element represents a coordinate of the observations\nf: value of the observations minus the background estimate (m-by-1 array).   (see note)\nlen: correlation length\nepsilon2: error variance of the observations (normalized by the error variance of the background field). epsilon2 can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a difference error variance and their errors are decorrelated) or a matrix (all observations can have a difference error variance and their errors can be correlated). If epsilon2 is a scalar, it is thus the inverse of the signal-to-noise ratio.\nnl: number of testing points around the current value of L. 1 means one additional point on both sides of the current L. 0 is allowed and means the parameter is not optimised.\nne: number of testing points around the current value of epsilon2. 0 is allowed as for nl\nmethod: cross validation estimator method 1: full CV  2: sampled CV 3: GCV 0: automatic choice between the three possible ones, default value\nOptional input arguments specified via keyword arguments are the same as for divand\n\nOutput:\n\nbestfactorl: best estimate of the multiplication factor to apply to len\nbestfactore: best estimate of the multiplication factor to apply to epsilon2\ncvvales : the cross validation values calculated\nfactors : the tested multiplication factors\ncvinter : interpolated cv values for final optimisation\nX2Data, Y2Data : coordinates of sampled cross validation in L,epsilon2 space . Normally only used for debugging or plotting\nXi2D, Yi2D : coordinates of interpolated estimator . Normally only used for debugging or plotting\n\nThe output bestfactorl and bestfactore represent multiplication factors which should be applied to L and epsilon2.\n\nThe len and epsilon2 provided should be close the real one as the tests will be performed around.\n\n\n\n"
},

{
    "location": "index.html#divand.empiriccovar",
    "page": "divand.jl documentation",
    "title": "divand.empiriccovar",
    "category": "function",
    "text": "distx,covar,corr,varx,count = empiriccovar(x,v,distbin,mincount;\n                          maxpoints = 10000,\n                          distfun = (xi,xj) -> sqrt(sum(abs2,xi-xj)))\n\nCompute the covariance, correlation and variance of a cloud of data points with the value v (a vector) and the location x (a tuple of vectors) grouped by distance. Random pairs are choosen and grouped by their distance (computed by distfun) in bins defined by distbin. The function try to fill at least mincount of data points in each bin but always stop after considering maxpoints pairs.\n\n\n\n"
},

{
    "location": "index.html#divand.fithorzlen",
    "page": "divand.jl documentation",
    "title": "divand.fithorzlen",
    "category": "function",
    "text": "lenz,dbinfo = divand.fithorzlen(x,value,z)\n\nDetermines the horizontal correlation length lenz based on the measurments value at the location x (tuple of 3 vectors corresponding to longitude, latitude and depth).\n\nOptional arguments:\n\ndistbin (default 0:0.5:5): distance bins to compute the empirical covariance functions\nmincount (default 50): minimum pairs per bins\nmaxpoints (default 10000): maximum number of pairs\nnmean (default 100): number of empirical covariance to average\nminlen (default 0.1): minimum correlation length-scale\nmaxlen (default 3): maximum correlation length-scale\nsmoothz (default 100): spatial filter for the correlation scale\nsearchz (default 300): vertical search distance\nmaxnsamp (default 5000): maximum number of samples\nlimitlen (default false): limit correlation length by mean distance between observations\n\n\n\n"
},

{
    "location": "index.html#divand.fitvertlen",
    "page": "divand.jl documentation",
    "title": "divand.fitvertlen",
    "category": "function",
    "text": "lenz,dbinfo = divand.fitvertlen(x,value,z)\n\nSee also divand.fithorzlen\n\n\n\n"
},

{
    "location": "index.html#divand.lengraddepth",
    "page": "divand.jl documentation",
    "title": "divand.lengraddepth",
    "category": "function",
    "text": "RL = lengraddepth(pmn,h, L;\n                  h2 = h,\n                  hmin = 0.001\n                  )\n\nCreate the relative correlation length-scale field RL based on the bathymetry  h and the metric pmn (tuple of arrays). Effectively the correlation-length  scale is close to zero if the relative bathymetry gradients (|∇h|/h) are smaller  than the length-scale L (in consistent units as pmn).\n\nR_L = 1 / (1 + L |∇h| / max(h2,hmin))\n\nPer default h2 is equal to h. The depth h must be positive. hmin must  have the same units as h (usually meters).\n\n\n\n"
},

{
    "location": "index.html#divand.divand_cvestimator",
    "page": "divand.jl documentation",
    "title": "divand.divand_cvestimator",
    "category": "function",
    "text": "theta = divand_cvestimator(s,residual)\n\nComputes the cross validation estimator (d-hatd)^T mathbf R^-1 (d-hatd)  ( mathbf 1^T mathbf R^-1 mathbf 1) where the hatd is the analysis not using a data point.\n\n\n\n"
},

{
    "location": "index.html#divand.weight_RtimesOne",
    "page": "divand.jl documentation",
    "title": "divand.weight_RtimesOne",
    "category": "function",
    "text": " weights = weight_RtimesOne(x,len)\n\nCompute the weight of the observations at location x to reduce the influence  of locally clustered data.  x is a tuple with n elements. Every element  represents a coordinate of the observations. len is a tuple of arrays representing the correlation length. len[i] is the correlation length in the  i-th dimension.\n\n\n\n"
},

{
    "location": "index.html#divand.Rtimesx!",
    "page": "divand.jl documentation",
    "title": "divand.Rtimesx!",
    "category": "function",
    "text": "Rtimesx!(coord,LS,x,Rx)\n\nGaussian type R matix in ndim dimensions applied to vector x of length  ndata. The Gaussian scale differs in each direction k : LS[k] Coordinates of point i are coord[i,1],coord[i,2],...,coord[i,ndim] To avoid an ndata² complexity a grid is set up first so as to allow only to calculate covarances when distances are smaller than 3*LS\n\nAdapted from DIVA3D/src/Fortran/Util/Rtimesx_weighting.f90\n\n\n\n"
},

{
    "location": "index.html#Parameter-optimization-1",
    "page": "divand.jl documentation",
    "title": "Parameter optimization",
    "category": "section",
    "text": "divand.fit_isotropic\ndivand.fit\ndivand.divand_cv\ndivand.empiriccovar\ndivand.fithorzlen\ndivand.fitvertlen\ndivand.lengraddepth\ndivand.divand_cvestimator\ndivand.weight_RtimesOne\ndivand.Rtimesx!"
},

{
    "location": "index.html#divand.Vocab.@urn_str",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.@urn_str",
    "category": "macro",
    "text": "urn\"SDN:x:y:z\'\n\nResolve a SeaDataNet URN (Uniform Resource Name) using https://www.seadatanet.org/urnurl/\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.CFVocab",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.CFVocab",
    "category": "type",
    "text": "collection = Vocab.CFVocab()\ncollection = Vocab.CFVocab(url = url)\n\nCreate a Dict-like object represeting the NetCDF CF Standard Name vocabulary. If the url is not provided then current CF Standard Name list http://cfconventions.org/Data/cf-standard-names/current/src/cf-standard-name-table.xml is used. Individual standard names are retirved by indexing which return an object of the type CFEntry:\n\ncollection = Vocab.CFVocab()\nentry = collection[\"sea_water_temperature\"]\n\n\n\n"
},

{
    "location": "index.html#Base.haskey-Tuple{divand.Vocab.CFVocab,Any}",
    "page": "divand.jl documentation",
    "title": "Base.haskey",
    "category": "method",
    "text": "bool = haskey(collection::CFVocab,stdname)\n\nReturn true if stdname is part of the NetCDF CF Standard Name vocabulary collection.\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.SDNCollection",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.SDNCollection",
    "category": "function",
    "text": "collection = SDNCollection(name)\n\nOpen the SeaDataNet collection with the name name at the URL http://www.seadatanet.org/urnurl/collection/ The collection can be indexed with brackets using the identifier.\n\nusing divand\ncollection = Vocab.SDNCollection(\"P01\")\nconcept = collection[\"PSALPR01\"]\n@show Vocab.prefLabel(concept)\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.prefLabel",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.prefLabel",
    "category": "function",
    "text": "s = Vocab.prefLabel(c::Vocab.Concept)\n\nReturn the preferred label of a concept c\n\n\n\ns = Vocab.prefLabel(urn::AbstractString)\n\nReturn the preferred label of a concept usings it URN (Uniform Resource Name)\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.altLabel",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.altLabel",
    "category": "function",
    "text": "s = Vocab.altLabel(c::Vocab.Concept)\n\nReturn the alternative label of a concept c\n\n\n\ns = Vocab.altLabel(urn::AbstractString)\n\nReturn the alternative label of a concept usings it URN (Uniform Resource Name)\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.notation",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.notation",
    "category": "function",
    "text": "s = Vocab.notation(c::Vocab.Concept)\n\nReturn the identifier of a concept c\n\n\n\ns = Vocab.notation(urn::AbstractString)\n\nReturn the identifier of a concept usings it URN (Uniform Resource Name)\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.definition",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.definition",
    "category": "function",
    "text": "s = Vocab.definition(c::Vocab.Concept)\n\nReturn the definition of a concept c\n\n\n\ns = Vocab.definition(urn::AbstractString)\n\nReturn the definition of a concept usings it URN (Uniform Resource Name)\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.resolve",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.resolve",
    "category": "function",
    "text": "entry = Vocab.resolve(urn)\n\nResolve a SeaDataNet URN (Uniform Resource Name) and returns the corresponding EDMO entry or Vocabulary concept. For example:\n\nconcept = Vocab.resolve(\"SDN:P021:current:TEMP\")\n\n\n\n"
},

{
    "location": "index.html#Base.find-Tuple{divand.Vocab.Concept,Any,Any}",
    "page": "divand.jl documentation",
    "title": "Base.find",
    "category": "method",
    "text": "find(c::Concept,name,collection)\n\nReturn a list of related concepts in the collection collection. name can be the string \"related\", \"narrower\", \"broader\".\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.description",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.description",
    "category": "function",
    "text": "str = description(entry::CFEntry)\nstr = canonical_units(entry::CFEntry)\n\nReturn the description or the canonical units of the `entry`.\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.canonical_units",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.canonical_units",
    "category": "function",
    "text": "str = description(entry::CFEntry)\nstr = canonical_units(entry::CFEntry)\n\nReturn the description or the canonical units of the `entry`.\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.splitURL",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.splitURL",
    "category": "function",
    "text": "collection,tag,key = Vocab.splitURL(url)\n\nSplit a concept URL into collection, tag and key. url must finishe with a slash.\n\n\n\n"
},

{
    "location": "index.html#Vocabulary-1",
    "page": "divand.jl documentation",
    "title": "Vocabulary",
    "category": "section",
    "text": "divand.Vocab.@urn_str\ndivand.Vocab.CFVocab\nBase.haskey(collection::divand.Vocab.CFVocab,stdname)\ndivand.Vocab.SDNCollection\ndivand.Vocab.prefLabel\ndivand.Vocab.altLabel\ndivand.Vocab.notation\ndivand.Vocab.definition\ndivand.Vocab.resolve\ndivand.Vocab.find(c::divand.Vocab.Concept,name,collection)\ndivand.Vocab.description\ndivand.Vocab.canonical_units\ndivand.Vocab.splitURL"
},

{
    "location": "index.html#Internal-API-or-advanced-usage-1",
    "page": "divand.jl documentation",
    "title": "Internal API or advanced usage",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#divand.statevector",
    "page": "divand.jl documentation",
    "title": "divand.statevector",
    "category": "type",
    "text": "Initialize structure for packing and unpacking given their mask.\n\nsv = statevector_init((mask1, mask2, ...))\n\nInitialize structure for packing and unpacking multiple variables given their corresponding land-sea mask.\n\nInput:   mask1, mask2,...: land-sea mask for variable 1,2,... Sea grid points correspond to one and land grid points to zero.     Every mask can have a different shape.\n\nOutput:   sv: structure to be used with statevector_pack and statevector_unpack.\n\nNote: see also statevector_pack, statevector_unpack\n\nAuthor: Alexander Barth, 2009,2017 <a.barth@ulg.ac.be> License: GPL 2 or later\n\n\n\n"
},

{
    "location": "index.html#divand.pack",
    "page": "divand.jl documentation",
    "title": "divand.pack",
    "category": "function",
    "text": "Pack a series of variables into a vector under the control of a mask.\n\nx = pack(sv,(var1, var2, ...))\n\nPack the different variables var1, var2, ... into the vector x where sv is a statevector. Only sea grid points are retained.\n\nInput:   sv: structure generated by statevector_init.   var1, var2,...: variables to pack (with the same shape as the corresponding masks).\n\nOutput:   x: vector of the packed elements. The size of this vector is the number of elements of all masks equal to 1.\n\nNotes: If var1, var2, ... have an additional trailing dimension, then this dimension is assumed to represent the different ensemble members. In this case x is a matrix and its last dimension is the number of ensemble members.\n\n\n\n"
},

{
    "location": "index.html#divand.unpack",
    "page": "divand.jl documentation",
    "title": "divand.unpack",
    "category": "function",
    "text": "Unpack a vector into different variables under the control of a mask.\n\nvar1, var2, ... = unpack(sv,x) var1, var2, ... = unpack(sv,x,fillvalue)\n\nUnpack the vector x into the different variables var1, var2, ... where sv is a statevector.\n\nInput:   sv: structure generated by statevector_init.   x: vector of the packed elements. The size of this vector is the number of elements equal to 1     in all masks.\n\nOptional input parameter:   fillvalue: The value to fill in var1, var2,... where the masks correspond to a land grid point. The default is zero.\n\nOutput:   var1, var2,...: unpacked variables.\n\nNotes: If x is a matrix, then the second dimension is assumed to represent the different ensemble members. In this case, var1, var2, ... have also an additional trailing dimension.\n\n\n\n"
},

{
    "location": "index.html#Base.sub2ind",
    "page": "divand.jl documentation",
    "title": "Base.sub2ind",
    "category": "function",
    "text": "ind = statevector_sub2ind(sv,subscripts)\n\nCompute from a tuple of subscripts the linear index in the packed state vector. The first element of the subscript indicates the variable index and the remaining the spatial subscripts.\n\n\n\n"
},

{
    "location": "index.html#Base.ind2sub",
    "page": "divand.jl documentation",
    "title": "Base.ind2sub",
    "category": "function",
    "text": "subscripts = ind2ind(sv,index)\n\nCompute from linear index in the packed state vector a tuple of subscripts. The first element of the subscript indicates the variable index and the remaining the spatial subscripts.\n\n\n\n"
},

{
    "location": "index.html#Base.length",
    "page": "divand.jl documentation",
    "title": "Base.length",
    "category": "function",
    "text": "number of points per node it is always zero for non-leaf nodes\n\n\n\n"
},

{
    "location": "index.html#State-vector-1",
    "page": "divand.jl documentation",
    "title": "State vector",
    "category": "section",
    "text": "divand.statevector\ndivand.pack\ndivand.unpack\nBase.sub2ind\nBase.ind2sub\nBase.length"
},

{
    "location": "index.html#Constraints-1",
    "page": "divand.jl documentation",
    "title": "Constraints",
    "category": "section",
    "text": "divand_constr_fluxes\ndivand_constr_constcoast"
},

{
    "location": "index.html#divand.ODVspreadsheet.listSDNparams",
    "page": "divand.jl documentation",
    "title": "divand.ODVspreadsheet.listSDNparams",
    "category": "function",
    "text": "p = listSDNparam(ODVData)\n\nReturn a list of SeaDataNet P01 parameters in a ODV spreadsheet `ODVData`.\n\n\n\n"
},

{
    "location": "index.html#divand.ODVspreadsheet.load",
    "page": "divand.jl documentation",
    "title": "divand.ODVspreadsheet.load",
    "category": "function",
    "text": " profiles,lons,lats,depths,times,ids = load(T,fnames,datanames;\n    qv_flags = [divand.ODVspreadsheet.GOOD_VALUE,\n                divand.ODVspreadsheet.PROBABLY_GOOD_VALUE],\n    nametype = :P01,\n    qvlocalname = \"QV:SEADATANET\")\n\nLoad all profiles in all file from the array fnames corresponding to one of the parameter names datanames. If nametype is :P01 (default), the datanames are P01 vocabulary names with the SDN prefix. If nametype is :localname, then they are the ODV column header without units. For example if the column header is Water body salinity [per mille], then datenames should be [\"Water body salinity\"]. The resulting vectors have the data type T (expect times and ids which are vectors of DateTime and String respectively). Only values matching the quality flag qv_flags are retained. qv_flags is a vector of Strings (based on http://vocab.nerc.ac.uk/collection/L20/current/, e.g. \"1\" means \"good value\"). One can also use the constants these constants (prefixed with divand.ODVspreadsheet.):\n\nqvlocalname is the column name to denote quality flags. It is assumed that the quality flags follow immediatly the data column.\n\nconstant value\nNO_QUALITY_CONTROL \"0\"\nGOOD_VALUE \"1\"\nPROBABLY_GOOD_VALUE \"2\"\nPROBABLY_BAD_VALUE \"3\"\nBAD_VALUE \"4\"\nCHANGED_VALUE \"5\"\nVALUE_BELOW_DETECTION \"6\"\nVALUE_IN_EXCESS \"7\"\nINTERPOLATED_VALUE \"8\"\nMISSING_VALUE \"9\"\nVALUE_PHENOMENON_UNCERTAIN \"A\"\n\nIf the ODV does not contain a semantic header (e.g. for the aggregated ODV files), then local names must be used.\n\njulia> data,lon,lat,depth,time,ids = divand.ODVspreadsheet.load(Float64,[\"data_from_med_profiles_non-restricted_v2.txt\"],\n      [\"Water body salinity\"]; nametype = :localname );\n\nNo checks are done if the units are consistent.\n\n\n\n profiles,lons,lats,depths,times,ids = load(T,dir,P01names)\n\nLoad all ODV files under the directory dir corresponding the one of the parameter names P01names. The resulting vectors have the data type T (expect times and ids which are vectors of DateTime and String respectively).\n\nNo checks are done if the units are consistent.\n\n\n\n"
},

{
    "location": "index.html#divand.ODVspreadsheet.localnames",
    "page": "divand.jl documentation",
    "title": "divand.ODVspreadsheet.localnames",
    "category": "function",
    "text": "list = localnames(sheet,P01name)\n\nReturn a list list of all local names mapping to the specified P01name in the ODV spreadsheet sheet without the prefix \"SDN:LOCAL:\".\n\n\n\nlist = localnames(sheet)\n\nReturn a list list of all local names  in the ODV spreadsheet sheet without the prefix \"SDN:LOCAL:\" in the order as they appear in the ODV file.\n\n\n\n"
},

{
    "location": "index.html#divand.ODVspreadsheet.Spreadsheet",
    "page": "divand.jl documentation",
    "title": "divand.ODVspreadsheet.Spreadsheet",
    "category": "type",
    "text": "Define composite type that will contain:\n\nthe metadata (dictionary),\nSDN parameter mapping (dictionary)\nthe column labels (array) and\nthe profiles (array of arrays).\n\n\n\n"
},

{
    "location": "index.html#divand.ODVspreadsheet.loadprofile",
    "page": "divand.jl documentation",
    "title": "divand.ODVspreadsheet.loadprofile",
    "category": "function",
    "text": " data,data_qv,obslon,obslat,obsdepth,obsdepth_qv,obstime,obstime_qv,EDMO,LOCAL_CDI_ID =\n loadprofile(T,sheet,iprofile,dataname; nametype = :P01)\n\nLoad a iprofile-th profile from the ODV spreadsheet sheet of the parameter dataname. If nametype is :P01 (default), the dataname is the P01 vocabulary name with the SDN prefix. If nametype is :localname, then it is the ODV column header.  The resulting vectors have the data type T (expect the quality flag and obstime) .\n\n\n\n"
},

{
    "location": "index.html#divand.ODVspreadsheet.loaddataqv",
    "page": "divand.jl documentation",
    "title": "divand.ODVspreadsheet.loaddataqv",
    "category": "function",
    "text": "data,data_qv = loaddataqv(sheet,profile,locname,fillvalue; fillmode = :repeat)\n\nThe same as loaddata, but now the quality flag are also loaded.\n\nprofile[i][j] is the j-th column of the i-th row of a profile. profile[i,j] is the i-th column of the j-th row of a profile.\n\n\n\n"
},

{
    "location": "index.html#divand.ODVspreadsheet.SDNparse!",
    "page": "divand.jl documentation",
    "title": "divand.ODVspreadsheet.SDNparse!",
    "category": "function",
    "text": "SDNparse!(col,fillmode,fillvalue,data)\n\nParse the list of String col into the corresponding data type of the vector data. Empty values are either replaced by fillvalue (if fillmode is :fill) or the previous value if repeated (if fillmode is :repeat)\n\n\n\n"
},

{
    "location": "index.html#divand.ODVspreadsheet.colnumber",
    "page": "divand.jl documentation",
    "title": "divand.ODVspreadsheet.colnumber",
    "category": "function",
    "text": "cn = colnumber(sheet,localname)\n\nReturn the column number cn of the first column with the local name localname (without the prefix \"SDN:LOCAL:\") in the ODV spreadsheet sheet.\n\n\n\n"
},

{
    "location": "index.html#divand.ODVspreadsheet.nprofiles",
    "page": "divand.jl documentation",
    "title": "divand.ODVspreadsheet.nprofiles",
    "category": "function",
    "text": "n = nprofiles(ODVData)\n\nReturn the number of profiles in a ODV Spreadsheet ODVData loaded by readODVspreadsheet.\n\n\n\n"
},

{
    "location": "index.html#ODV-files-1",
    "page": "divand.jl documentation",
    "title": "ODV files",
    "category": "section",
    "text": "divand.ODVspreadsheet.listSDNparams\ndivand.ODVspreadsheet.load\ndivand.ODVspreadsheet.localnames\ndivand.ODVspreadsheet.Spreadsheet\ndivand.ODVspreadsheet.loadprofile\ndivand.ODVspreadsheet.loaddataqv\ndivand.ODVspreadsheet.SDNparse!\ndivand.ODVspreadsheet.colnumber\ndivand.ODVspreadsheet.nprofiles"
},

{
    "location": "index.html#divand.sparse_interp",
    "page": "divand.jl documentation",
    "title": "divand.sparse_interp",
    "category": "function",
    "text": "H,out = sparse_interp(mask,I)\n\nCreate interpolation matrix from mask and fractional indexes I.\n\nInput:   mask: 0 invalid and 1 valid points (n-dimensional array)   I: fractional indexes (2-dim array n by mi, where mi is the number of points to interpolate) Ouput:   H: sparse matrix with interpolation coefficients   out: true if value outside of grid   outbbox: 1 if outise bouding box   onland: 1 if point touches land (where mask == 0)\n\n\n\n"
},

{
    "location": "index.html#divand.sparse_interp_g",
    "page": "divand.jl documentation",
    "title": "divand.sparse_interp_g",
    "category": "function",
    "text": "sparse_interp(x,mask,xi) Interpolate from x onto xi\n\n\n\n"
},

{
    "location": "index.html#divand.sparse_gradient",
    "page": "divand.jl documentation",
    "title": "divand.sparse_gradient",
    "category": "function",
    "text": "Sparse operator for a gradient. Dx1,Dx2,...,Dxn = sparse_gradient(mask,pmn) Form the gradient using finite differences in all n-dimensions Input:   mask: binary mask delimiting the domain. 1 is inside and 0 outside.         For oceanographic application, this is the land-sea mask.   pmn: scale factor of the grid. Output:   Dx1,Dx2,...,Dxn: operators represeting a gradient along     different dimensions\n\n\n\n"
},

{
    "location": "index.html#divand.sparse_diff",
    "page": "divand.jl documentation",
    "title": "divand.sparse_diff",
    "category": "function",
    "text": "Sparse operator for differentiation.\n\ndiffx = sparse_diff(sz1,m,cyclic)\n\nSparse operator for differentiation along dimension m for \"collapsed\" matrix of the size sz1.\n\nInput:   sz1: size of rhs   m: dimension to differentiate   cyclic: true if domain is cyclic along dimension m. False is the   default value\n\n\n\n"
},

{
    "location": "index.html#divand.matfun_trim",
    "page": "divand.jl documentation",
    "title": "divand.matfun_trim",
    "category": "function",
    "text": "T = matfun_trim(sz1,m)\n\nCreate an operator which trim first and last row (or column) in The field is a \"collapsed\" matrix of the size sz1. m is the dimension  to trim.\n\n\n\n"
},

{
    "location": "index.html#divand.matfun_stagger",
    "page": "divand.jl documentation",
    "title": "divand.matfun_stagger",
    "category": "function",
    "text": "S = matfun_stagger(sz1,m,cyclic)\n\nCreate an operator for staggering a field in dimension m. The field is a \"collapsed\" matrix of the size sz1.\n\nInput:   sz1: size of rhs   m: dimension to stagger   cyclic: true if domain is cyclic along dimension m. False is the   default value\n\n\n\n"
},

{
    "location": "index.html#divand.matfun_diff",
    "page": "divand.jl documentation",
    "title": "divand.matfun_diff",
    "category": "function",
    "text": "Operator for differentiation.\n\ndiffx = matfun_diff(sz1,m,cyclic)\n\nOperator for differentiation along dimension m for \"collapsed\" matrix of the size sz1.\n\nInput:   sz1: size of rhs   m: dimension to differentiate   cyclic: true if domain is cyclic along dimension m. False is the   default value\n\n\n\n"
},

{
    "location": "index.html#divand.matfun_shift",
    "page": "divand.jl documentation",
    "title": "divand.matfun_shift",
    "category": "function",
    "text": "Operator shifting a field in a given dimension.\n\nfunction S = matfun_shift(sz1,m,cyclic)\n\nOperator shifting a field in the dimension m. The field is a \"collapsed\" matrix of the size sz1.\n\nInput:   sz1: size of rhs   m: dimension to shift   cyclic: true if domain is cyclic along dimension m. False is the     default value\n\n\n\n"
},

{
    "location": "index.html#Operators-1",
    "page": "divand.jl documentation",
    "title": "Operators",
    "category": "section",
    "text": "divand.sparse_interp\ndivand.sparse_interp_g\ndivand.sparse_gradient\ndivand.sparse_diff\ndivand.matfun_trim\ndivand.matfun_stagger\ndivand.matfun_diff\ndivand.matfun_shift"
},

{
    "location": "index.html#divand.Quadtrees.QT",
    "page": "divand.jl documentation",
    "title": "divand.Quadtrees.QT",
    "category": "type",
    "text": "quadtree (of the higher-dimensional equivalent) T the type of the coordinates TA the type of the attributes N number of dimensions\n\n\n\n"
},

{
    "location": "index.html#divand.Quadtrees.rsplit!",
    "page": "divand.jl documentation",
    "title": "divand.Quadtrees.rsplit!",
    "category": "function",
    "text": "recursive split\n\n\n\n"
},

{
    "location": "index.html#divand.Quadtrees.add!",
    "page": "divand.jl documentation",
    "title": "divand.Quadtrees.add!",
    "category": "function",
    "text": "sucess = add!(qt,x,attrib,max_cap = 10) Add point x with the attribute attrib to the quadtree qt. sucess is true if xis within the bounds of the quadtree node qt (otherwise  false and the point has not been added)\n\n\n\n"
},

{
    "location": "index.html#divand.Quadtrees.within",
    "page": "divand.jl documentation",
    "title": "divand.Quadtrees.within",
    "category": "function",
    "text": "attribs = within(qt,min,max)\n\nSearch all points within a bounding box defined by the vectors min and max.\n\n\n\n"
},

{
    "location": "index.html#divand.Quadtrees.bitget",
    "page": "divand.jl documentation",
    "title": "divand.Quadtrees.bitget",
    "category": "function",
    "text": "Test if the n-th bit in a is set. The least significant bit is n = 1.\n\n\n\n"
},

{
    "location": "index.html#divand.Quadtrees.inside",
    "page": "divand.jl documentation",
    "title": "divand.Quadtrees.inside",
    "category": "function",
    "text": "         x1\n\n+–––––+    |          |   |   +      |   |   y      |   +–––––+  x0\n\n\n\n"
},

{
    "location": "index.html#divand.Quadtrees.intersect",
    "page": "divand.jl documentation",
    "title": "divand.Quadtrees.intersect",
    "category": "function",
    "text": "Test of the rectanges defined by x0,x1  and y0,y1 intersects              x1   +–––––+    |          |   |   +–––––+ y1   |   |      |   |   +–––––+   |  x0   |          |       |          |       +–––––+      y0\n\n\n\n"
},

{
    "location": "index.html#divand.Quadtrees.split!",
    "page": "divand.jl documentation",
    "title": "divand.Quadtrees.split!",
    "category": "function",
    "text": "split a single node\n\n\n\n"
},

{
    "location": "index.html#Quadtree-1",
    "page": "divand.jl documentation",
    "title": "Quadtree",
    "category": "section",
    "text": "divand.Quadtrees.QT\ndivand.Quadtrees.rsplit!\ndivand.Quadtrees.add!\ndivand.Quadtrees.within\ndivand.Quadtrees.bitget\ndivand.Quadtrees.inside\ndivand.Quadtrees.intersect\ndivand.Quadtrees.split!"
},

{
    "location": "index.html#divand.conjugategradient",
    "page": "divand.jl documentation",
    "title": "divand.conjugategradient",
    "category": "function",
    "text": "x,success,niter = conjugategradient(fun!,b)\n\nSolve a linear system with the preconditioned conjugated-gradient method: A x = b where A is a symmetric positive defined matrix and b is a vector.  Equivalently the solution x minimizes the cost function  J(x) = ½ xᵀ A x - bᵀ x.\n\nThe function fun!(x,fx) computes fx which is equal to  A*x. For example:\n\nfunction fun!(x,fx)\n    fx[:] = A*x\nend\n\nNote that the following code will NOT work, because a new array fx would be created and it would not be passed back to the caller.\n\nfunction fun!(x,fx)\n    fx = A*x # bug!\nend\n\nThe function fun! works in-place to reduce the amount of memory allocations.\n\nOptional input arguments\n\nx0: starting vector for the interations\ntol: tolerance on  |Ax-b| / |b|\nmaxit: maximum of interations\npc!: the preconditioner. The functions pc(x,fx) computes fx = M⁻¹ x (the inverse of M times x) where M is a symmetric positive defined matrix. Effectively, the system E⁻¹ A (E⁻¹)ᵀ (E x) = E⁻¹ b is solved for (E x) where E Eᵀ = M. Ideally, M should this be similar to A, so that E⁻¹ A (E⁻¹)ᵀ is close to the identity matrix. The function pc! should be implemented in a similar way than fun! (see above).\n\nOutput\n\nx: the solution\nsuccess: true if the interation converged (otherwise false)\nniter: the number of iterations\n\n\n\n"
},

{
    "location": "index.html#divand.pc_none!",
    "page": "divand.jl documentation",
    "title": "divand.pc_none!",
    "category": "function",
    "text": "pc_none!(x,fx)\n\nDummy call-back function when no preconditioner is used. fx will be equal to x.\n\n\n\n"
},

{
    "location": "index.html#divand.checksym",
    "page": "divand.jl documentation",
    "title": "divand.checksym",
    "category": "function",
    "text": "xAy, yATx = checksym(n,fun!)\n\nCheck if the the function fun! represents a symmetric matrix when applied on  random vectors of size n. \n\n\n\n"
},

{
    "location": "index.html#Conjugate-gradient-1",
    "page": "divand.jl documentation",
    "title": "Conjugate gradient",
    "category": "section",
    "text": "divand.conjugategradient\ndivand.pc_none!\ndivand.checksym"
},

{
    "location": "index.html#divand.divand_laplacian",
    "page": "divand.jl documentation",
    "title": "divand.divand_laplacian",
    "category": "function",
    "text": "Create the laplacian operator.\n\nLap = divand_laplacian(mask,pmn,nu,iscyclic)\n\nForm a Laplacian using finite differences  assumes that gradient is zero at \"coastline\"\n\nInput:    mask: binary mask delimiting the domain. 1 is inside and 0 outside.          For oceanographic application, this is the land-sea mask.    pmn: scale factor of the grid.    nu: diffusion coefficient of the Laplacian       field of the size mask or cell arrays of fields\n\nOutput:    Lap: sparce matrix represeting a Laplacian\n\n\n\n"
},

{
    "location": "index.html#divand.divand_obscovar",
    "page": "divand.jl documentation",
    "title": "divand.divand_obscovar",
    "category": "function",
    "text": "R = divand_obscovar(epsilon2,m)\n\nCreate a matrix representing the observation error covariance R of size m x m.\n\nIf epsilon2 is a scalar, then R = epsilon2 * I If epsilon2 is a vector, then R = diag(epsilon2) If epsilon2 is a matrix, then R = epsilon2\n\n\n\n"
},

{
    "location": "index.html#divand.divand_adaptedeps2",
    "page": "divand.jl documentation",
    "title": "divand.divand_adaptedeps2",
    "category": "function",
    "text": "factor = divand_adaptedeps2(s,fi);\n\nInput:\n\ns: structure returned by divandrun\nfi: analysis returned by divandrun\n\nOutput:\n\nfactor : multiplicative factor to apply to epsilon2\n\nUsing Deroziers adaptive approach provides a multiplicative factor for the current epsilon2 value so that factor*epsilon2 is a better estimate of the R matrix. If you cannot use divandrun but use divandgo, the latter provides automatically this pamater as result.\n\n\n\n"
},

{
    "location": "index.html#divand.divand_diagHKobs",
    "page": "divand.jl documentation",
    "title": "divand.divand_diagHKobs",
    "category": "function",
    "text": "Computes the diagonal terms of the so called hat-matrix HK, using the already solved analysis and it structure s. Warning: might take some time\n\nThis version only uses the real data (not those related to additional constraints)\n\ndiagonalterms = divand_diagHKobs(s);\n\n\n\n"
},

{
    "location": "index.html#divand.divand_residual",
    "page": "divand.jl documentation",
    "title": "divand.divand_residual",
    "category": "function",
    "text": "dataresidual = divand_residual(s,fi)\n\nComputes the generalized residual yo - H xa  using the analysis on the grid  fi and the solution structure s.\n\n\n\n"
},

{
    "location": "index.html#divand.divand_addc",
    "page": "divand.jl documentation",
    "title": "divand.divand_addc",
    "category": "function",
    "text": "s = divand_addc(s,c)\n\nAdd a constraint c to the cost function defined by s. The structure s is typically created by divand_background and the contrain c  has the following fields: R (a covariance matrix), H (extraction operator) and  yo (specified value for the constrain). The added contrain Jc(x) is quadratic and has the following structure.\n\nJc(x) = (H x - yo)ᵀ R⁻¹ (H x - yo)\n\n\n\n"
},

{
    "location": "index.html#divand.divand_erroratdatapoints",
    "page": "divand.jl documentation",
    "title": "divand.divand_erroratdatapoints",
    "category": "function",
    "text": "Computes the error at the real data locations using the analysis structure s\n\nerrorvariance = divand_erroratdatapoints(s);\n\n\n\n"
},

{
    "location": "index.html#divand.divand_iBpHtiRHx!",
    "page": "divand.jl documentation",
    "title": "divand.divand_iBpHtiRHx!",
    "category": "function",
    "text": "\n\n"
},

{
    "location": "index.html#divand.divand_GCVKii",
    "page": "divand.jl documentation",
    "title": "divand.divand_GCVKii",
    "category": "function",
    "text": "Computes an estimate of the mean value of the diagonal of HK using GCV and the already solved analysisand it structure s\n\nKii = divand_GCVKii(s);\n\n\n\n"
},

{
    "location": "index.html#divand.divand_fittocpu",
    "page": "divand.jl documentation",
    "title": "divand.divand_fittocpu",
    "category": "function",
    "text": "stepsize,overlapping,isdirect = divand_fittocpu(Lpmnrange,gridsize,latercsteps,moddim=[]);\n\nCreates a list of windows for subsequent domain decomposition\n\nAlso calculates already the subsampling steps csteps for the preconditionners\n\nInput:\n\nLpmnrange:\ngridsize: number of points in each direction (size(mask))\nmoddim:\n\nOutput:\n\n\n\n"
},

{
    "location": "index.html#divand.divand_background",
    "page": "divand.jl documentation",
    "title": "divand.divand_background",
    "category": "function",
    "text": "Form the inverse of the background error covariance matrix. s = divand_background(mask,pmn,Labs,alpha,moddim) Form the inverse of the background error covariance matrix with finite-difference operators on a curvilinear grid\n\nInput:\n\nmask: binary mask delimiting the domain. 1 is inside and 0 outside.       For oceanographic application, this is the land-sea mask.\npmn: scale factor of the grid.\nLabs: correlation length\nalpha: a dimensional coefficients for norm, gradient, laplacian,...    alpha is usually [1,2,1] in 2 dimensions.\n\nOutput:\n\ns: stucture containing\ns.iB: inverse of the background error covariance\ns.L: spatial average correlation length\ns.n: number of dimenions\ns.coeff: scaling coefficient such that the background variance diag(inv(iB)) is one far away from the boundary.\n\n\n\n"
},

{
    "location": "index.html#divand.divand_obs",
    "page": "divand.jl documentation",
    "title": "divand.divand_obs",
    "category": "function",
    "text": "s = divand_obs(s,xi,x,R,I)\n\nInclude the constrain from the observations. It is assumed that the each coordinate depends only on one index. If this is not the case, then matrix I must be provided.\n\nInput:   s: structure created by divand_background   xi: coordinates of observations (tuple of vectors)   x: coordinates of grid (tuple of arrays)   R: obs. error covariance matrix (normalized)   I (optional): fractional indexes of location of observation     within the grid\n\nOutput:   s: structure to be used by divand_factorize\n\nNote make sure not to mix Float32 and Float64 for divand_constrain.\n\n\n\n"
},

{
    "location": "index.html#divand.divand_bc_stretch",
    "page": "divand.jl documentation",
    "title": "divand.divand_bc_stretch",
    "category": "function",
    "text": "\n\n"
},

{
    "location": "index.html#divand.divand_diagHK",
    "page": "divand.jl documentation",
    "title": "divand.divand_diagHK",
    "category": "function",
    "text": "Computes the diagonal terms of the so called hat-matrix HK, using the already solved analysis and it structure s. Warning: might take some time\n\ndiagonalterms = divand_diagHK(s);\n\n\n\n"
},

{
    "location": "index.html#divand.divand_kernel",
    "page": "divand.jl documentation",
    "title": "divand.divand_kernel",
    "category": "function",
    "text": "mu,K,len_scale = divand_kernel(n,alpha)\n\nReturn the analytical kernel and normalization factor.\n\nAnalytical (normalized) kernels K for infinite domain in dimension n and for coefficients alpha and normalization factor mu.\n\nK(r) is the kernel function (function of the normalized distance r), len_scale is the distance at which K(len_scale) = 0.6019072301972346 (which is besselk(1,1))\n\n\n\n"
},

{
    "location": "index.html#divand.divand_residualobs",
    "page": "divand.jl documentation",
    "title": "divand.divand_residualobs",
    "category": "function",
    "text": "Computes the residual yo- H xa  only at real data points using the analysis on the grid fi and the solution structure s\n\ndataresidual = divand_residualobs(s,fi);\n\n\n\n"
},

{
    "location": "index.html#divand.divand_aexerr",
    "page": "divand.jl documentation",
    "title": "divand.divand_aexerr",
    "category": "function",
    "text": "aexerr,Bref,fa,sa = divand_aexerr(mask,pmn,xi,x,f,len,epsilon2;...);\n\nInput: same as for divandrun\n\nmask: binary mask delimiting the domain. true is inside and false outside. For oceanographic application, this is the land-sea mask.\npmn: scale factor of the grid. pmn is a tuple with n elements. Every      element represents the scale factor of the corresponding dimension. Its      inverse is the local resolution of the grid in a particular dimension.\nxi: tuple with n elements. Every element represents a coordinate of the final grid on which the observations are interpolated\nx: tuple with n elements. Every element represents a coordinate of the observations\nf: value of the observations minus the background estimate (m-by-1 array).   (see note)\nlen: correlation length\nepsilon2: error variance of the observations (normalized by the error variance of the background field). epsilon2 can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a difference error variance and their errors are decorrelated) or a matrix (all observations can have a difference error variance and their errors can be correlated). If epsilon2 is a scalar, it is thus the inverse of the signal-to-noise ratio.\n\nOptional input arguments specified as keyword arguments also as for divand\n\nOutput:\n\naexerr: the almost exact error\nBref: the background error for error scaling by background aexerr./Bref\nfa: the analysis (with low impact fake data): DO NOT USE UNLESS YOU KNOW WHAT YOU ARE DOING\nsa: the associated structure\n\nCompute a variational analysis of arbitrarily located observations to calculate the almost exact error \n\n\n\n"
},

{
    "location": "index.html#divand.divand_cpme",
    "page": "divand.jl documentation",
    "title": "divand.divand_cpme",
    "category": "function",
    "text": "cpme = divand_cpme(mask,pmn,xi,x,f,len,epsilon2;...);\n\nInput: Same as for divandrun\n\nmask: binary mask delimiting the domain. true is inside and false outside. For oceanographic application, this is the land-sea mask.\npmn: scale factor of the grid. pmn is a tuple with n elements. Every      element represents the scale factor of the corresponding dimension. Its      inverse is the local resolution of the grid in a particular dimension.\nxi: tuple with n elements. Every element represents a coordinate of the final grid on which the observations are interpolated\nx: tuple with n elements. Every element represents a coordinate of the observations\nf: value of the observations minus the background estimate (m-by-1 array).   (see note)\nlen: correlation length\nepsilon2: error variance of the observations (normalized by the error variance of the background field). epsilon2 can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a difference error variance and their errors are decorrelated) or a matrix (all observations can have a difference error variance and their errors can be correlated). If epsilon2 is a scalar, it is thus the inverse of the signal-to-noise ratio.\nkeywords : undocumented for the moment how to use iterative solver with coarser grid as preconditionner. see divandjog for csteps, lmask and alphapcparameters\n\nOptional input arguments specified as keyword arguments also as for divand\n\nOutput:\n\ncpme: the clever poor mans error\n\nPerform an n-dimensional variational analysis of the observations f located at the coordinates x. The array cpme represent the error field at the grid defined by the coordinates xi and the scales factors pmn. If you cannot run divandrun you can use divandgo with error field calculation :cpme\n\n\n\n"
},

{
    "location": "index.html#divand.divand_cpme_go",
    "page": "divand.jl documentation",
    "title": "divand.divand_cpme_go",
    "category": "function",
    "text": "erri = divand_cpme_go(mask,pmn,xi,x,f,len,epsilon2; ...);\n\nInput:\n\nSame arguments as divandrun with in addition\nMEMTOFIT=: keyword controlling how to cut the domain depending on the memory remaining available for inversion (not total memory)\nRTIMESONESCALES= : if you provide a tuple of length scales, data are weighted differently depending on the numbers of neighbours they have. See weight_RtimesOne for details \n\nOutput:\n\nerri: relative error field using the clever poor man\'s error approach. Result on the same grid as fi. `\n\nONLY USE THIS VERSION IF YOU CANNOT RUN divandgo with :cmpe activated (or directly divand_cpme if you can run divandrun)\n\n\n\n"
},

{
    "location": "index.html#divand.divand_datainboundingbox",
    "page": "divand.jl documentation",
    "title": "divand.divand_datainboundingbox",
    "category": "function",
    "text": "divand_datainboundingbox(xi,x,f;Rmatrix=())\n\nInput:\n\nxi: tuple with n elements. Every element represents a coordinate   of the final grid on which the observations are interpolated\n\nx: tuple with n elements. Every element represents a coordinate of the observations\nf: value of the observations\nRmatrix: error variance of the observations (normalized by the error variance of the background field). epsilon2 can be a scalar (all observations have the same error variance and their errors are decorrelated), a vector (all observations can have a difference error variance and their errors are decorrelated) or a matrix (all observations can have a difference error variance and their errors can be correlated). If epsilon2 is a scalar, it is thus the inverse of the signal-to-noise ratio.\n\nOutput:\n\nxn: tuple with n elements. Every element represents a coordinate of   the observations which falls in the bounding box defined by xi fn: the corresponding data indexes: the indexes in the original array retained  Rn: the new error variance \n\n\n\n"
},

{
    "location": "index.html#divand.divand_Lpmnrange",
    "page": "divand.jl documentation",
    "title": "divand.divand_Lpmnrange",
    "category": "function",
    "text": "Lpmnrange = divand_Lpmnrange(pmn,len);\n\nIn each direction searches for the minimum and maximum value of the length scale times the metric in this diretion\n\nSi it basically looks at the worst and the best resolution found in the grid\n\nInput:\n\npmn: scale factor of the grid. pmn is a tuple with n elements. Every      element represents the scale factor of the corresponding dimension. Its      inverse is the local resolution of the grid in a particular dimension.\nlen: correlation length\n\nOutput:\n\nLpmnrange: Array of range tuples (minimum and maximum of L times metric)\n\n\n\n"
},

{
    "location": "index.html#divand.divand_pc_sqrtiB",
    "page": "divand.jl documentation",
    "title": "divand.divand_pc_sqrtiB",
    "category": "function",
    "text": "Compute a preconditioner using the Cholesky decomposition.\n\n[M1,M2] = divand_pc_michol(iB,H,R)\n\nCompute preconditioner matrices M1 and M2 based on the Cholesky decomposition of iB. The matrices H and R are not used. M2 is the transpose of M1 for this preconditioner.\n\n\n\n"
},

{
    "location": "index.html#divand.divand_pc_none",
    "page": "divand.jl documentation",
    "title": "divand.divand_pc_none",
    "category": "function",
    "text": "fun = divand_pc_none(iB,H,R)\n\nDummy function for requiring that no preconditioner is used in divand.\n\nSee also: diavnd_pc_sqrtiB\n\n\n\n"
},

{
    "location": "index.html#divand.divand_GCVKiiobs",
    "page": "divand.jl documentation",
    "title": "divand.divand_GCVKiiobs",
    "category": "function",
    "text": "Computes an estimate of the mean value of the diagonal of HK using GCV and the already solved analysis and it structure s\n\nOnly using real data locations\n\nKii = divand_GCVKiiobs(s);\n\n\n\n"
},

{
    "location": "index.html#divand.divand_cutter",
    "page": "divand.jl documentation",
    "title": "divand.divand_cutter",
    "category": "function",
    "text": "windowlist,csteps,lmask,alphapc = divand_cutter(Lpmnrange,gridsize,moddim,MEMTOFIT);\n\nCreates a list of windows for subsequent domain decomposition\n\nAlso calculates already the subsampling steps csteps for the preconditionners\n\nas well as the mask lmask to apply to the length scales in the preconditionner, allowing to reduce\n\nthe problem size\n\nInput:\n\nLpmnrange:\ngridsize: number of points in each direction (size(mask))\nmoddim:\n\nOutput:\n\nwindowlist: Array of tuples (iw1 iw2 ...)\ncsteps : Array of steps for the coarse grid preconditionner\nlmask : Array of multiplication factors for length scale of preconditionner\nalphapc : Norm defining coefficients for preconditionner\n\n\n\n"
},

{
    "location": "index.html#divand.divand_qc",
    "page": "divand.jl documentation",
    "title": "divand.divand_qc",
    "category": "function",
    "text": "qcvalues = divand_qc(fi,s,method);\n\nInput:\n\nfi : interpolated field from a divandrun execution\ns: corresponding structure returned by divand\nmethod : optional argument, which describes the method to be used:\n\n1  as for standard cross validation,  3  as for GCV,  4  with CV estimator to be used outside the routine,  5  Poor man\'s GCV using data instead of random vector,  0  automatic selection of method.\n\nOutput\n\nqcvalues: quality check values, one for each data point.\n\nThe higher the value, the more suspect a data point is. Absolute values of qcvalues might be not robust when analysis parameters are uncertain. The ranking is however quite robust.\n\nIf you cannot run divandrun but use divandgo (which does not provide a structure s at the output), the latter provides qcvalues if you call divandgo with a keyword parameter QCMETHOD=\n\n\n\n"
},

{
    "location": "index.html#divand.divand_solve!",
    "page": "divand.jl documentation",
    "title": "divand.divand_solve!",
    "category": "function",
    "text": "Solve the variational problem.\n\nfi = divand_solve(s)\n\nDerive the analysis based on all contraints included in s and using the observations yo\n\nInput:   s: structure created by divand_factorize   fi0: starting point for iterative primal methods   f0: starting point for the iterative dual method\n\nbtrunc: the value at which the stored value of s.iB was truncated and needs to be completed on the fly using jmBix\n\nOutput:   fi: analyzed field\n\n\n\n"
},

{
    "location": "index.html#divand.divand_sampler",
    "page": "divand.jl documentation",
    "title": "divand.divand_sampler",
    "category": "function",
    "text": "samplesteps = divand_sampler(pmn,len);\n\nDefines steps for sub-sampling in the discrete grid which would still allow to resolve the provided lengthscales\n\nInput:\n\npmn: scale factor of the grid. pmn is a tuple with n elements. Every      element represents the scale factor of the corresponding dimension. Its      inverse is the local resolution of the grid in a particular dimension.\nlen: correlation length\n\nOutput:\n\nsamplesteps: vector of integers with steps in subsampling [1 2 4 1] means every grid point in x direction, every fifth in y etc\n\n\n\n"
},

{
    "location": "index.html#divand.divandjog",
    "page": "divand.jl documentation",
    "title": "divand.divandjog",
    "category": "function",
    "text": "Compute a variational analysis of arbitrarily located observations.\n\nfi,s = divandjog(mask,pmn,xi,x,f,len,epsilon2,csteps,lmask; alphapc=[1,2,1], otherargs...);\n\nPerform an n-dimensional variational analysis of the observations f located at the coordinates x. The array fi represent the interpolated field at the grid defined by the coordinates xi and the scales factors pmn.\n\nInput:\n\nSame parameters as for divarun.       * Two additional parameters:               * csteps: array of ndims values defining the sampling steps for the preconditionner               * lmask: array of ndims mutilplications factors for length scales       * One additional optiional parameter               * alphapc: The coefficients for the norm used in the preconditionner\n\nOutput:\n\nfi: the analysed field\ns: structure with an array s.P representing the analysed error covariance\n\nNote:\n\n\n\n"
},

{
    "location": "index.html#divand.divand_background_components",
    "page": "divand.jl documentation",
    "title": "divand.divand_background_components",
    "category": "function",
    "text": "Form the different components of the background error covariance matrix.\n\niB = divand_background_components(s,D,alpha; kwargs...)\n\nCompute the components of the background error covariance matrix iB_ and their sum based on alpha (the a-dimensional coefficients for norm, gradient, laplacian,...).\n\nIf the optional arguments contains btrunc, the calculation of iB is limited to the term up and including alpha[btrunc]\n\n\n\n"
},

{
    "location": "index.html#divand.stats",
    "page": "divand.jl documentation",
    "title": "divand.stats",
    "category": "function",
    "text": "meanx,stdx = stats(sumx,sumx2,N)\n\nComputes the mean meanx and the standard deviation stdx from the sum (sumx) and the sum of squares (sumx2) from N numbers.\n\n\n\nmeanx,meany,stdx,stdy,covar,corr = stats(sumx,sumx2,sumy,sumy2,sumxy,N)\n\nComputes the mean meanx and the standard deviation stdx from the sum (sumx) and the sum of squares (sumx2) from N numbers and similarly for the variable y. The function computes also the Pearson correlation corr and covariance covar between x and y.\n\n\n\n"
},

{
    "location": "index.html#divand.statpos",
    "page": "divand.jl documentation",
    "title": "divand.statpos",
    "category": "function",
    "text": "ulon,ulat,meanval,stdval,count = statpos(val,lon,lat)\n\nReturn unique positions (ulon, ulat) as well their mean, standard deviation and count of the vector of observations val located at the positions lon and lat.\n\n\n\n"
},

{
    "location": "index.html#divand.blkdiag",
    "page": "divand.jl documentation",
    "title": "divand.blkdiag",
    "category": "function",
    "text": "concatenate diagonal matrices\n\n\n\n"
},

{
    "location": "index.html#Base.findfirst",
    "page": "divand.jl documentation",
    "title": "Base.findfirst",
    "category": "function",
    "text": "findfirst(c::Concept,name,collection)\n\nReturn the first related concepts in the collection collection. name can be the string \"related\", \"narrower\", \"broader\".\n\n\n\n"
},

{
    "location": "index.html#divand.formatsize",
    "page": "divand.jl documentation",
    "title": "divand.formatsize",
    "category": "function",
    "text": "display size as a string \n\n\n\n"
},

{
    "location": "index.html#divand.interp!",
    "page": "divand.jl documentation",
    "title": "divand.interp!",
    "category": "function",
    "text": "interp!(xi,fi,x,f)\n\nInterpolate field fi (n-dimensional array) defined at xi (tuble of n-dimensional arrays or vectors) onto grid x (tuble of n-dimensional arrays). The interpolated field is stored in f. The grid in xi must be align with the axis (e.g. produced by divand.ndgrid).\n\n\n\n"
},

{
    "location": "index.html#divand.ufill",
    "page": "divand.jl documentation",
    "title": "divand.ufill",
    "category": "function",
    "text": "cfilled = ufill(c,valex)\n\nReplace values in c equal to valex by averages of surrounding points.\n\n\n\nufill(c::Array{T,2},mask::AbstractArray{Bool}) where T\n\nmask is true where c is valid.\n\n\n\n"
},

{
    "location": "index.html#divand.jmBix",
    "page": "divand.jl documentation",
    "title": "divand.jmBix",
    "category": "function",
    "text": "\n\n"
},

{
    "location": "index.html#divand.cgradient",
    "page": "divand.jl documentation",
    "title": "divand.cgradient",
    "category": "function",
    "text": "hx,hy = cgradient(pmn,h)\n\n\n\n"
},

{
    "location": "index.html#divand.fzero",
    "page": "divand.jl documentation",
    "title": "divand.fzero",
    "category": "function",
    "text": "fzero(f,x0,x1,eps; maxiter = Inf) find the zero of the function f between x0 and x1 assuming x0 < x1 at a precision eps.\n\n\n\n"
},

{
    "location": "index.html#divand.localize_separable_grid",
    "page": "divand.jl documentation",
    "title": "divand.localize_separable_grid",
    "category": "function",
    "text": "Derive fractional indices on a separable grid.\n\nI = localize_separable_grid(xi,mask,x)\n\nxi and x are a tuples, e.g. x1,x2 = ndgrid(2 * collect(1:5),collect(1:6)) x = (x1,x2)\n\nDerive fractional indices where xi are the points to localize in the separable grid x (every dimension in independent on other dimension). The output I is an n-by-m array where n number of dimensions and m number of observations. The correspond element of I is negative if xi is outside of the grid defined by x.\n\n\n\n"
},

{
    "location": "index.html#divand.decompB!",
    "page": "divand.jl documentation",
    "title": "divand.decompB!",
    "category": "function",
    "text": "work1, work2: size of mask\n\nSymmetric matrix\n\nSB = √(β) (1 + α L)^(nmax/2) W^{-1}\n\nwhere W is the volumne of the corresponding grid cell. The background error covariance matrix B is SB W SB\n\n\n\n"
},

{
    "location": "index.html#divand.varanalysis",
    "page": "divand.jl documentation",
    "title": "divand.varanalysis",
    "category": "function",
    "text": "Variational analysis similar to 3D-var\n\nInput:\n\nx0: start vector for iteration, at output it is the last state of the     iteration. Note that x0 is related to the analysis xa by       xa = SB^1/2 * W^1/2 * xa\n\n| x + W^1/2 * SB^1/2 * H\' * (R  (H * SB^1/2 * W^1/2 * x ))   -   W^1/2 SB^{1/2} * H\' * (R  yo) |       <     tol * s.sv.n / length(yo)  * | W^1/2 SB^{1/2} * H\' * (R  yo) |\n\nKernel is the solution of the n-dimensional diffusion equation\n\n∂c/∂t =  ∇ ⋅ (D ∇ c)\n\nn-dimensional Green’s function\n\nG(x,x\',t) = (4πDt)^(-n/2)  exp( - |x -x\'|² / (4Dt))\n\nG(x,x\',t) = det(D)^(-1/2) (4π t)^(-n/2)  exp( - (x -x\')ᵀ D⁻¹ (x -x\')ᵀ / (4t))\n\nhttp://www.rpgroup.caltech.edu/~natsirt/aph162/diffusion_old.pdf\n\n\n\n"
},

{
    "location": "index.html#divand.len_harmonize",
    "page": "divand.jl documentation",
    "title": "divand.len_harmonize",
    "category": "function",
    "text": "Len = len_harmonise(len,mask) Produce a tuple of arrays of the correlation length len which can be either a scalar (homogeneous and isotropic case), a tuple of scalar (homogeneous case) or already a tuple of arrays (general case). The the later case the size of the arrays are veryfied.\n\n\n\n"
},

{
    "location": "index.html#divand.alpha_default",
    "page": "divand.jl documentation",
    "title": "divand.alpha_default",
    "category": "function",
    "text": "neff, alpha = alpha_default(Labs,alpha)\n\nReturn a default value of alpha.\n\n\n\n"
},

{
    "location": "index.html#divand.ncfile",
    "page": "divand.jl documentation",
    "title": "divand.ncfile",
    "category": "function",
    "text": "divand_save(ds,filename,xyi,fi,varname;\n                  ncvarattrib = Dict(), ncglobalattrib = Dict(), ...)\n\nSave the result of the analysis in a NetCDF file .\n\nInput arguments\n\nds: the NetCDF dataset \nfilename: the name of the NetCDF file\nmask: binary mask delimiting the domain. true is inside and false outside. For oceanographic application, this is the land-sea mask where sea is true and land is false.\nxyi: tuple with n elements. Every element represents a coordinate of the final grid on which the observations are interpolated\nfi: the analysed field\nvarname: the name of the NetCDF variable\n\nOptional arguments:\n\nncglobalattrib: a dictionary with the global attributes\nncvarattrib: a dictionary with the variable attributes\nrelerr: relative error\ntimeorigin: time origin for the time units attribute (default is 1900-01-01 00:00:00)\n\n\n\n"
},

{
    "location": "index.html#divand.writeslice",
    "page": "divand.jl documentation",
    "title": "divand.writeslice",
    "category": "function",
    "text": "ncvar, ncvar_relerr, ncvar_Lx, fi, relerr, index)\n\nWhite a slice of data in a NetCDF given by the index index. The variable relerr can be nothing.\n\n\n\n"
},

{
    "location": "index.html#divand.encodeWMSStyle",
    "page": "divand.jl documentation",
    "title": "divand.encodeWMSStyle",
    "category": "function",
    "text": "encode parameters as key-value separated by : and +\n\n\n\n"
},

{
    "location": "index.html#divand.loadoriginators",
    "page": "divand.jl documentation",
    "title": "divand.loadoriginators",
    "category": "function",
    "text": "db = loadoriginators(fname)\n\nLoad the CDI list from the file fname (zip with a csv file, or csv file directly).\n\n\n\n"
},

{
    "location": "index.html#Utility-functions-1",
    "page": "divand.jl documentation",
    "title": "Utility functions",
    "category": "section",
    "text": "divand.divand_laplacian\ndivand.divand_obscovar\ndivand.divand_adaptedeps2\ndivand.divand_diagHKobs\ndivand.divand_residual\ndivand.divand_addc\ndivand.divand_erroratdatapoints\ndivand.divand_iBpHtiRHx!\ndivand.divand_GCVKii\ndivand.divand_fittocpu\ndivand.divand_background\ndivand.divand_obs\ndivand.divand_bc_stretch\ndivand.divand_diagHK\ndivand.divand_kernel\ndivand.divand_residualobs\ndivand.divand_aexerr\ndivand.divand_cpme\ndivand.divand_cpme_go\ndivand.divand_datainboundingbox\ndivand.divand_Lpmnrange\ndivand.divand_pc_sqrtiB\ndivand.divand_pc_none\ndivand.divand_GCVKiiobs\ndivand.divand_cutter\ndivand.divand_qc\ndivand.divand_solve!\ndivand.divand_sampler\ndivand.divandjog\ndivand.divand_background_components\ndivand.stats\ndivand.statpos\ndivand.blkdiag\nBase.findfirst\ndivand.formatsize\ndivand.interp!\ndivand.ufill\ndivand.jmBix\ndivand.cgradient\ndivand.fzero\ndivand.localize_separable_grid\ndivand.decompB!\ndivand.varanalysis\ndivand.len_harmonize\ndivand.alpha_default\ndivand.ncfile\ndivand.writeslice\ndivand.encodeWMSStyle\ndivand.loadoriginators"
},

{
    "location": "index.html#Examples-1",
    "page": "divand.jl documentation",
    "title": "Examples",
    "category": "section",
    "text": "To run the example, you need to install PyPlot. In the folder examples of divand, you can run e.g. the example divand_simple_example_1D.jl by issuing:# cd(\"/path/to/divand/examples\")\ninclude(\"divand_simple_example_1D.jl\")Replace /path/to/divand/ by the installation directory of divand which is the output of Pkg.dir(\"divand\") if you installed divand using Julias package manager."
},

{
    "location": "index.html#Performance-considerations-1",
    "page": "divand.jl documentation",
    "title": "Performance considerations",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Tuning-the-domain-decomposition-1",
    "page": "divand.jl documentation",
    "title": "Tuning the domain decomposition",
    "category": "section",
    "text": "The functions diva3d and divandgo split the domain into overlapping subdomains to reduce the required amount of memory. In some circumstances (in particular few vertical levels), this can unnecessarily degrade the performance. The CPU time of the analysis can be improved by increasing the diva3d option memtofit from 3 (default) to higher values (as long as one does not run out of memory). If this parameter is set to a very high value then the domain decomposition is effectively disabled."
},

{
    "location": "index.html#Multiple-CPU-system-1",
    "page": "divand.jl documentation",
    "title": "Multiple CPU system",
    "category": "section",
    "text": "Per default julia tries to use all CPUs on your system when doing matrix operations. The number of CPUs is controlled by the call to BLAS.set_num_threads. Using multiple CPUs can result in overhead and it can be beneficial to reduce the number of CPUs:BLAS.set_num_threads(2)"
},

{
    "location": "index.html#Information-for-developers-1",
    "page": "divand.jl documentation",
    "title": "Information for developers",
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

{
    "location": "index.html#Troubleshooting-1",
    "page": "divand.jl documentation",
    "title": "Troubleshooting",
    "category": "section",
    "text": "If the installation of a package fails, it is recommended to update the local copy of the package list by issuing Pkg.update() to make sure that Julia knows about the latest version of these packages and then to re-try the installation of the problematic package. Julia calls the local copy of the packge list METADATA. For example to retry the installation of EzXML issue the following command:Pkg.update()\nPkg.add(\"EzXML\")"
},

{
    "location": "index.html#No-plotting-window-appears-1",
    "page": "divand.jl documentation",
    "title": "No plotting window appears",
    "category": "section",
    "text": "If the following command doesn\'t produce any figureusing PyPlot\nplot(1:10)A possible solution is to modify the backend: this is done by editing the python configuration file matplotlibrc. The location of this file is obtained in python with:import matplotlib\nmatplotlib.matplotlib_fnameUnder Linux, this returns \'~/.config/matplotlib/matplotlibrc\'. To use the TkAgg backend, add the following to the file:backend      : TkAggThe matplotlibrc need to be created if it does not exists."
},

{
    "location": "index.html#Julia-cannot-connect-to-GitHub-on-Windows-7-and-Windows-Server-2012-1",
    "page": "divand.jl documentation",
    "title": "Julia cannot connect to GitHub on Windows 7 and Windows Server 2012",
    "category": "section",
    "text": "Cloning METADATA or downloading a julia packages fails with:GitError(Code:ECERTIFICATE, Class:OS, , user cancelled certificate checks: )The problem is that Windows 7 and Windows Server 2012 uses outdated encryption protocols. The solution is to run the \"Easy fix\" tool from the Microsoft support page"
},

{
    "location": "index.html#MbedTLS.jl-does-not-install-on-Windows-7-1",
    "page": "divand.jl documentation",
    "title": "MbedTLS.jl does not install on Windows 7",
    "category": "section",
    "text": "The installion of MbedTLS.jl fails with the error message:INFO: Building MbedTLS                                                                                                                                    \nInfo: Downloading https://github.com/quinnj/MbedTLSBuilder/releases/download/v0.6/MbedTLS.x86_64-w64-mingw32.tar.gz to C:\\Users\\Jeremy\\.julia\\v0.6\\MbedTLS\n\\deps\\usr\\downloads\\MbedTLS.x86_64-w64-mingw32.tar.gz...                                                                                                  \nException setting \"SecurityProtocol\": \"Cannot convert null to type \"System.Net.SecurityProtocolType\" due to invalid enumeration values. Specify one of th\ne following enumeration values and try again. The possible enumeration values are \"Ssl3, Tls\".\"                                                           \nAt line:1 char:35                                                                                                                                         \n+ [System.Net.ServicePointManager]:: <<<< SecurityProtocol =                                                                                              \n    + CategoryInfo          : InvalidOperation: (:) [], RuntimeException                                                                                  \n    + FullyQualifiedErrorId : PropertyAssignmentException                                                                                                 \n    [...]See also the issue https://github.com/JuliaWeb/MbedTLS.jl/issues/133The solution is to install the Windows Management Framework 4.0."
},

{
    "location": "index.html#EzXML.jl-cannot-be-installed-on-RedHat-6-1",
    "page": "divand.jl documentation",
    "title": "EzXML.jl cannot be installed on RedHat 6",
    "category": "section",
    "text": "The zlib library of RedHat 6, is slightly older than the library which EzXML.jl and libxml2 requires.To verify this issue, you can type in JuliaLibdl.dlopen(joinpath(Pkg.dir(\"EzXML\"),\"deps/usr/lib/libxml2.so\"))It should not return an error message. On Redhat 6.6, the following error message is returned:ERROR: could not load library \"/home/username/.julia/v0.6/EzXML/deps/usr/lib/libxml2.so\"\n\n/lib64/libz.so.1: version `ZLIB_1.2.3.3\' not found (required by /home/divahs1/.julia/v0.6/EzXML/deps/usr/lib/libxml2.so)\n\nStacktrace:\n\n [1] dlopen(::String, ::UInt32) at ./libdl.jl:97 (repeats 2 times)However, the following command should work: LD_LIBRARY_PATH=\"$HOME/.julia/v0.6/EzXML/deps/usr/lib/:$LD_LIBRARY_PATH\" julia --eval  \'print(Libdl.dlopen(joinpath(Pkg.dir(\"EzXML\"),\"deps/usr/lib/libxml2.so\"))\'Lukily, EzZML.jl includes a newer version of the zlib library, but it does not load the library automatically. (see also https://github.com/JuliaLang/julia/issues/7004 and https://github.com/JuliaIO/HDF5.jl/issues/97)To make Julia use this library, a user on RedHat 6 should always start Julia with:LD_LIBRARY_PATH=\"$HOME/.julia/v0.6/EzXML/deps/usr/lib/:$LD_LIBRARY_PATH\" juliaOne can also create script with the following content:#!/bin/bash\nexport LD_LIBRARY_PATH=\"$HOME/.julia/v0.6/EzXML/deps/usr/lib/:$LD_LIBRARY_PATH\"\nexec /path/to/bin/julia \"$@\"by replacing /path/to/bin/julia to the full path of your installation directory. The script should be marked executable and it can be included in your Linux search PATH environement variable. Julia can then be started by calling directly this script."
},

{
    "location": "index.html#The-DIVAnd-test-suite-fails-with-automatic-download-failed-1",
    "page": "divand.jl documentation",
    "title": "The DIVAnd test suite fails with automatic download failed",
    "category": "section",
    "text": "Running Pkg.test(\"divand\") fails with the error:automatic download failed (error: 2147500036)The test suite will download some sample data. You need to have internet access and run the test function from a directory with write access.You can change the directory to your home directory with the julia command cd(homedir()).You can check the current working directory with:pwd()"
},

{
    "location": "index.html#METADATA-cannot-be-updated-1",
    "page": "divand.jl documentation",
    "title": "METADATA cannot be updated",
    "category": "section",
    "text": "Pkg.update fails with the error message METADATA cannot be updated.If you have git installed, you can issue the command:cd ~/.julia/v0.6/METADATA\ngit reset --hardand then in Julia run Pkg.update() again.If this does not work, then, you can also delete ~/.julia (https://github.com/JuliaLang/julia/issues/18651#issuecomment-347579521) and in Julia enter Pkg.init(); Pkg.update()."
},

{
    "location": "index.html#Convert-error-in-divand_obs-1",
    "page": "divand.jl documentation",
    "title": "Convert error in divand_obs",
    "category": "section",
    "text": "The full error message:MethodError: Cannot `convert` an object of type divand.divand_constrain{Float32,Diagonal{Float64},SparseMatrixCSC{Float64,Int64}} to an object of type divand.divand_constrain{Float64,TR,TH} where TH<:(AbstractArray{#s370,2} where #s370<:Number) where TR<:(AbstractArray{#s371,2} where #s371<:Number)\nThis may have arisen from a call to the constructor divand.divand_constrain{Float64,TR,TH} where TH<:(AbstractArray{#s370,2} where #s370<:Number) where TR<:(AbstractArray{#s371,2} where #s371<:Number)(...),\nsince type constructors fall back to convert methods.The solution is to use the same type of all input parameters: all Float32 or all Float64."
},

{
    "location": "index.html#Monthlist-issue-1",
    "page": "divand.jl documentation",
    "title": "Monthlist issue",
    "category": "section",
    "text": "Using comments inside list can lead to unexpected results.This monthlist = [\n       [1,2,3]\n       #[4,5,6]\n       ]should be written as monthlist = [\n       [1,2,3]\n       ]"
},

{
    "location": "index.html#Error-in-the-factorisation-1",
    "page": "divand.jl documentation",
    "title": "Error in the factorisation",
    "category": "section",
    "text": "The following messageBase.LinAlg.PosDefException(95650)followed by the stack-trace starting with: julia Stacktrace:  [1] #cholfact!#8(::Float64, ::Function, ::Base.SparseArrays.CHOLMOD.Factor{Float64}, ::Base.SparseArrays.CHOLMOD.Sparse{Float64}) at ./sparse/cholmod.jl:1360  .................  [9] divandrun(::BitArray{3}, ::Tuple{Array{Float64,3},Array{Float64,3},Array{Float64,3}}, ::Tuple{Array{Float64,3},Array{Float64,3},Array{Float64,3}}, ::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}}, ::Array{Float64,1}, ::Tuple{Array{Float64,3},Array{Float64,3},Array{Float64,3}}, ::Float64) at /home/ctroupin/.julia/v0.6/divand/src/divandrun.jl:147might be due to a wrong choice in the analysis parameters, for example a too long  correlation length."
},

]}
