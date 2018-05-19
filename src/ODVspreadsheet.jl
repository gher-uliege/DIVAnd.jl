module ODVspreadsheet

#using Logging
using StringEncodings

# Set logging level(DEBUG, INFO, WARNING, ERROR or CRITICAL)
#loglevel = WARNING
#Logging.configure(level=loglevel);
debug(msg) = nothing

# SeaDataNet Quality Flags
# http://vocab.nerc.ac.uk/collection/L20/current/

const NO_QUALITY_CONTROL = "0"
const GOOD_VALUE = "1"
const PROBABLY_GOOD_VALUE = "2"
const PROBABLY_BAD_VALUE = "3"
const BAD_VALUE = "4"
const CHANGED_VALUE = "5"
const VALUE_BELOW_DETECTION = "6"
const VALUE_IN_EXCESS = "7"
const INTERPOLATED_VALUE = "8"
const MISSING_VALUE = "9"
const VALUE_PHENOMENON_UNCERTAIN = "A"

"""
Define composite type that will contain:
* the metadata (dictionary),
* SDN parameter mapping (dictionary)
* the column labels (array) and
* the profiles (array of arrays).
"""

global Spreadsheet
type Spreadsheet
    metadata::Dict{String,String}
    # local name to a tuple of object and unit
    # //<subject>SDN:LOCAL:Chronological Julian Date</subject><object>SDN:P01::CJDY1101</object><units>SDN:P06::UTAA</units>
    SDN_parameter_mapping::Dict{String,Dict{String,String}}
    columnLabels::Array{SubString{String},1}
    # profileList[k][i,j]
    # is the k-th profile, the i-th column (parameter) and
    # the j-th row (often depth or time)
    profileList::Array{Any,1}
end

"""
    n = nprofiles(ODVData)

Return the number of profiles in a ODV Spreadsheet `ODVData` loaded by
`readODVspreadsheet`.
"""

nprofiles(ODVData) = length(ODVData.profileList)


"""
    The function will return a composite type that will store:
    1. The general metadata of the spreadsheet
    2. The labels of the columns
    3. The individual profiles

    Input

    *`datafile`: the path to an ODV spreadsheet file.
               The Path can be relative or absolute.

    Output

    *`ODVdata`: a "Spreadsheet" composite type.

"""
function readODVspreadsheet(datafile)

    # metadata will be stored in a dictionary
    # ODV doc: Comment lines start with two slashes  // as first two characters
    metadata = Dict{String, String}()
    SDN_parameter_mapping = Dict{String,Dict{String,String}}()

    # using the encoding Latin-1 degrates significantly the performance
    info("Reading data from file $(datafile)")
    open(datafile, "r") do f
        line = readline(f)

        # Byte Order Mark (BOM)
        if startswith(line,"\uFEFF")
            # ignore BOM
            line = line[2:end]
        end

        # Read the metadata (lines starting with //)
        # unfortunately, there are also some files with empty line in the
        # header

        while startswith(line,"//") || (length(line) == 0)
            # Identify metadata fields using regex
            # (name of the field is between < > and </ >)
            m = match(r"<(\w+)>(.+)</(\w+)>", line)

            if m != nothing
                debug("Match found")
                debug(m[1] * ": " * m[2])
                # Add key - value in the dictionnary
                metadata[String(m[1])] = String(m[2])
            end

            if line == "//SDN_parameter_mapping"
                line = readline(f);
                info("Length line SDN_parameter_mapping: $(length(line))")
                # The semantic descriptions are terminated by an empty comment
                # record (i.e. a record containing the // characters and nothing else)

                while line != "//"
                    @assert startswith(line,"//")

                    # split at < or >
                    parts = split(line[3:end],r"[<|>]",keep=false)
                    tmp = Dict(k => v for (k,v) in zip(parts[1:3:end],parts[2:3:end]))

                    subject = tmp["subject"]
                    delete!(tmp,"subject")
                    SDN_parameter_mapping[subject] = tmp
                    line = readline(f);
                end
            end
            line = readline(f);
        end

        # Read the column labels and set number of columns
        #ODV doc: must be the first non-comment line in the file
        #ODV doc: must provide columns for all mandatory meta-variables
        columnline = line
        columnLabels = split(chomp(columnline), '\t')
        debug("Column labels: $(columnLabels)");
        ncols = length(columnLabels);
        info("Total no. of columns (before selection): $ncols")

        # Discard columns that won't be used (should be extended)
        column2discard = ["QF", "Instrument Info",
                          "Reference", "Data set name", "Discipline",
                          "Category", "Variables measured",
                          "Data format", "Data format version",
                          "Data size", "Data set creation date",
                          "Measuring area type", "Track resolution",
                          "Track resolution unit", "Frequency",
                          "Frequency unit", "Cruise name",
                          "Alternative cruise name", "Cruise start date",
                          "Station name", "Alternative station name",
                          "Station start date"];

        goodcols = setdiff(columnLabels, column2discard);
        # Get indices of the good columns
        index2keep = findin(columnLabels, goodcols);
        ncols2 = length(index2keep);
        info("No. of columns after selection: $ncols2")

        # number of total lines
        pos = position(f)
        totallines = 0
        # count the lines "cleverly", without the comments
        # or the empty lines
        for row in eachline(f)

            if startswith(row,"//")
                # ignore lines starting with e.g.
                # //<History> ...
                continue
            end

            line = split(row, "\t");

            # some files have only white space on the last line
            if all.(isempty(line))
                continue
            end
            totallines += 1
        end
        info("Total no. of lines: $totallines")
        seek(f,pos)

        # load all data
        alldata = Array{SubString{String},2}(ncols2,totallines);
        ##info("Size (in GB) of data matrix (OLD): " * string(sizeof(alldataold)))
        i = 0

        for row in eachline(f; chomp = true)
            i = i+1

            if startswith(row,"//")
                # ignore lines starting with e.g.
                # //<History> ...
                continue
            end

            line = split(row, "\t")[index2keep];

            # some files have only white space on the last line
            if all.(isempty(line))
                continue
            end

            if length(line) != ncols2
                error("Expecting $(ncols2) columns but $(length(line)) found in line $(line) (line number $i).")
            end

            alldata[:,i] = line
        end

        info("Size of data matrix: " * string(size(alldata)))
        # info("Size (in GB) of data matrix: " * string(sizeof(alldata / 1024 / 1024)))
        # trim unused lines (comments, ...)
        #alldata = alldata[:,1:i]

        # count profiles
        nprofiles = 0
        for i = 1:size(alldata,2)
            # If the first value (Station) is not empty,
            # then it's a profile
            if alldata[1,i] != ""
                nprofiles += 1
            end
        end

        # count each profiles length
        profile_start_index = Vector{Int}(nprofiles+1)
        j = 0
        for i = 1:size(alldata,2)
            # If the first value (Station) is not empty,
            # then it's a profile
            if alldata[1,i] != ""
                j += 1
                profile_start_index[j] = i
            end
        end
        profile_start_index[end] = size(alldata,2)+1

        # split into profiles
        iprofile = 1
        profileList = Vector{Array{SubString{String},2}}(nprofiles)
        j = 1
        for i = 1:nprofiles
            profileList[i] = alldata[:,profile_start_index[i]:profile_start_index[i+1]-1]
        end

        info("No. of profiles in the file: " * string(nprofiles))
        ODVdata = Spreadsheet(metadata, SDN_parameter_mapping, columnLabels[index2keep], profileList)
        return ODVdata
    end
end

"""
    p = listSDNparam(ODVData)

    Return a list of SeaDataNet P01 parameters in a ODV spreadsheet `ODVData`.
"""
function listSDNparams(ODVData)
    return [d["object"] for (k,d) in ODVData.SDN_parameter_mapping]
end


function countP01(sheets)

    count_P01 = Dict{String,Int}();
    for sheet in sheets
        for p in listSDNparams(sheet)
            count_P01[p] = get(count_P01,p,0) + 1;
        end
    end

    stat = sort([(k,v) for (k,v) in count_P01]; by = kv -> kv[2])
    return stat
end

"""
    list = localnames(sheet,P01name)

Return a list `list` of all local names mapping to the specified `P01name` in the
ODV spreadsheet `sheet` without the prefix "SDN:LOCAL:".

"""
localnames(sheet,P01name) = String[replace(v,r"^SDN:LOCAL:","") for (v,d) in sheet.SDN_parameter_mapping if d["object"] == P01name]

"""
    list = localnames(sheet)

Return a list `list` of all local names  in the
ODV spreadsheet `sheet` without the prefix "SDN:LOCAL:" in the order as they
appear in the ODV file.
"""
localnames(sheet) = String[strip(split(l,'[')[1]) for l in sheet.columnLabels]

"""
    cn = colnumber(sheet,localname)

Return the column number `cn` of the first column with the local name
`localname` (without the prefix "SDN:LOCAL:") in the ODV spreadsheet `sheet`.
"""

colnumber(sheet,localname) = findfirst(localnames(sheet) .== localname)


"""
    dt = parsejd(t)

Convert a Chronological Julian Day Number to a DateTime object. The
reference value is taken from
[Chronological Julian Date](https://web.archive.org/web/20171129142108/https://www.hermetic.ch/cal_stud/chron_jdate.htm)

From the SDN standard:
"A real number representing the Chronological Julian Date, which is defined as the time
elapsed in days from 00:00 on January 1 st 4713 BC. ... "

The time origin is _not_ noon (12:00) on Monday, January 1, 4713 BC as for the Julia Date Number.
"""

parsejd(t) = DateTime(2007,2,10) + Dates.Millisecond(round(Int64,(t - 2454142.) * (24*60*60*1000)))


"""
    v = myparse(T,s)

Parse the string `s` as a type `T`. Unlike Julia's parse function
an error message contains the string `s` (which could not be parsed) for
debugging.
"""
function myparse(T,s)
    try
        return parse(T,s)
    catch err
        if isa(err,ArgumentError)
            throw(ArgumentError("unable to parse $(s) as $(T)"))
        else
            rethrow(err)
        end
    end
end


"""
    SDNparse!(col,fillmode,fillvalue,data)

Parse the list of String `col` into the corresponding data type
of the vector `data`. Empty values are either replaced by `fillvalue`
(if fillmode is :fill) or the previous value if repeated (if fillmode
is :repeat)
"""

function SDNparse!(col,fillmode,fillvalue,data)
    for i = 1:length(col)
        if col[i] == ""
            if (i > 1) && (fillmode == :repeat)

                # section 2.3
                # As metadata values are constant throughout a row_group it is
                # usual practice just to populate the first row.

                data[i] = data[i-1]
            else
                data[i] = fillvalue
            end
        else
            if eltype(data) <: AbstractString
                data[i] = col[i]
            else
                data[i] = myparse(eltype(data),col[i])
                #data[i] = parse(eltype(data),col[i])
            end
        end
    end

    return data
end


"""
    data = loaddata(sheet,profile,locname,fillvalue; fillmode = :repeat)

Load a single column referred by the local name `locname` in the profile
`profile` from the ODV spreadsheet `sheet`. Empty values are either replaced
by `fillvalue` (if fillmode is :fill) or the previous value if repeated (if fillmode
is :repeat)
"""

function loaddata(sheet,profile,locname,fillvalue::T; fillmode = :repeat) where T
    lenprof = size(profile,2)

    cn_data = colnumber(sheet,locname)
    data = Vector{T}(lenprof)

    if cn_data == 0
        data[:] = fillvalue
    else
        SDNparse!(profile[cn_data,:],fillmode,fillvalue,data)
    end

    return data
end

"""
    data,data_qv = loaddataqv(sheet,profile,locname,fillvalue; fillmode = :repeat)

The same as `loaddata`, but now the quality flag are also loaded.

profile[i][j] is the j-th column of the i-th row of a profile.
profile[i,j] is the i-th column of the j-th row of a profile.
"""

function loaddataqv(sheet,profile,locname,fillvalue::T;
                    fillmode = :repeat,
                    qvlocalname = "QV:SEADATANET"
                    ) where T
    locnames = localnames(sheet)
    lenprof = size(profile,2)

    cn_data = colnumber(sheet,locname)
    data = Vector{T}(lenprof)
    data_qv = fill("",(lenprof,))


    if cn_data == 0
        data[:] = fillvalue
        data_qv[:] = MISSING_VALUE
    else
        if cn_data > length(profile)
            @show profile,cn_data
            error("incomplete profile record")
        end

        SDNparse!(profile[cn_data,:],fillmode,fillvalue,data)

        if (cn_data < length(locnames)) && (locnames[cn_data+1] == qvlocalname)
            data_qv[:] = profile[cn_data+1,:]
        end
    end

    return data,data_qv
end

"""
     data,data_qv,obslon,obslat,obsdepth,obsdepth_qv,obstime,obstime_qv,EDMO,LOCAL_CDI_ID =
     loadprofile(T,sheet,iprofile,dataname; nametype = :P01)

Load a `iprofile`-th profile from the ODV spreadsheet `sheet` of the
parameter `dataname`. If `nametype` is `:P01` (default), the
dataname is the P01 vocabulary name with the SDN prefix. If nametype is
`:localname`, then it is the ODV column header.
 The resulting vectors have the data type `T`
(expect the quality flag and `obstime`) .
"""

function loadprofile(T,sheet,iprofile,dataname;
                     nametype = :P01,
                     qvlocalname = "QV:SEADATANET"
                     )
    const fillvalue = T(NaN)
    const filldate_jd = 0.
    const filldate = parsejd(filldate_jd)

    profile = sheet.profileList[iprofile]
    locnames = localnames(sheet)
    P01names = listSDNparams(sheet)

    localname =
        if nametype == :P01
            localnames(sheet,dataname)[1]
        elseif nametype == :localname
            dataname
        else
            error("nametype should be :P01 or :localname and not $(nametype)")
        end

    data,data_qv = loaddataqv(sheet,profile,localname,fillvalue;
                              fillmode = :fill,
                              qvlocalname = qvlocalname)
    sz = size(data)

    #cruise = loaddata(sheet,profile,"Cruise","")
    #station = loaddata(sheet,profile,"Station","")
    #ptype = loaddata(sheet,profile,"Type","")

    # look for EDMO_code and EDMO_CODE
    # is both are absent return an empty value

    EDMO =
        if "EDMO_code" in locnames
            loaddata(sheet,profile,"EDMO_code","")
        else
            loaddata(sheet,profile,"EDMO_CODE","")
        end

    LOCAL_CDI_ID = loaddata(sheet,profile,"LOCAL_CDI_ID","")

    lon = loaddata(sheet,profile,"Longitude",fillvalue)
    lat = loaddata(sheet,profile,"Latitude",fillvalue)

    depth = fill(fillvalue,sz)
    depth_qv = fill("",sz)
    if "SDN:P01::ADEPZZ01" in P01names
        localname_depth = localnames(sheet,"SDN:P01::ADEPZZ01")
        depth[:],depth_qv[:] = loaddataqv(sheet,profile,localname_depth,fillvalue;
                                          qvlocalname = qvlocalname)
    elseif "Depth" in locnames
        depth[:],depth_qv[:] = loaddataqv(sheet,profile,"Depth",fillvalue;
                                          qvlocalname = qvlocalname)
        # if "Depth reference" in locnames
        #     depthref = loaddata(sheet,profile,"Depth reference","unknown")


        #     unexpected_depthref = ((depthref .!= "mean sea level") .&
        #                            (depthref .!= "sea level"))

        #     if any(unexpected_depthref)
        #         @show depthref[unexpected_depthref]
        #     end
        # end
    end

    time = Vector{DateTime}(sz)
    time_qv = Vector{String}(sz)

    # chronological julian day
    if "SDN:P01::CJDY1101" in P01names
        locname_time = localnames(sheet,"SDN:P01::CJDY1101")[1]
        timedata,time_qv = loaddataqv(sheet,profile,locname_time,filldate_jd;
                                      qvlocalname = qvlocalname)
        time = parsejd.(timedata)
    elseif "SDN:P01::DTUT8601" in P01names
        # ISO8601 format, e.g. yyyy-mm-ddThh:mm:ss.sss

        locname_time = localnames(sheet,"SDN:P01::DTUT8601")
        time,time_qv = loaddataqv(sheet,profile,locname_time,filldate;
                                  qvlocalname = qvlocalname)
    else
        # hopefully not necessary
        for header in ["yyyy-mm-ddThh:mm:ss.sss",
                       "yyyy-mm-ddThh:mm:ss",
                       "yyyy-mm-ddThh:mm",
                       "yyyy-mm-ddThh",
                       "yyyy-mm-dd"]
            if header in locnames
                time,time_qv = loaddataqv(sheet,profile,header,filldate;
                                          qvlocalname = qvlocalname)
                break
            end
        end
    end

    return data,data_qv,lon,lat,depth,depth_qv,time,time_qv,EDMO,LOCAL_CDI_ID
end


function goodflag(obstime_qv,qv_flags)
    good_time = falses(size(obstime_qv))
    for flag in qv_flags
        good_time[:] =  good_time .| (obstime_qv .== flag)
    end

    # time quality flag can also be absent
    good_time[:] =  good_time .| (obstime_qv .== "")

    return good_time
end


"""
     profiles,lons,lats,depths,times,ids = load(T,fnames,datanames;
        qv_flags = [divand.ODVspreadsheet.GOOD_VALUE,
                    divand.ODVspreadsheet.PROBABLY_GOOD_VALUE],
        nametype = :P01,
        qvlocalname = "QV:SEADATANET")

Load all profiles in all file from the array `fnames` corresponding to
one of the parameter names `datanames`. If `nametype` is `:P01` (default), the
datanames are P01 vocabulary names with the SDN prefix. If nametype is
`:localname`, then they are the ODV column header without units. For example
if the column header is `Water body salinity [per mille]`, then `datenames`
should be `["Water body salinity"]`.
The resulting vectors have the data type `T` (expect `times` and `ids` which
are vectors of `DateTime` and `String` respectively). Only values matching the
quality flag `qv_flags` are retained. `qv_flags` is a vector of Strings
(based on http://vocab.nerc.ac.uk/collection/L20/current/, e.g. "1" means "good value").
One can also use the constants these constants (prefixed with
`divand.ODVspreadsheet.`):

`qvlocalname` is the column name to denote quality flags. It is assumed that the
quality flags follow immediatly the data column.

|                    constant | value |
|-----------------------------|-------|
|          NO_QUALITY_CONTROL |   "0" |
|                  GOOD_VALUE |   "1" |
|         PROBABLY_GOOD_VALUE |   "2" |
|          PROBABLY_BAD_VALUE |   "3" |
|                   BAD_VALUE |   "4" |
|               CHANGED_VALUE |   "5" |
|       VALUE_BELOW_DETECTION |   "6" |
|             VALUE_IN_EXCESS |   "7" |
|          INTERPOLATED_VALUE |   "8" |
|               MISSING_VALUE |   "9" |
|  VALUE_PHENOMENON_UNCERTAIN |   "A" |


If the ODV does not contain a semantic header (e.g. for the aggregated ODV files),
then local names must be used.

```julia-repl
julia> data,lon,lat,depth,time,ids = divand.ODVspreadsheet.load(Float64,["data_from_med_profiles_non-restricted_v2.txt"],
      ["Water body salinity"]; nametype = :localname );
```

No checks are done if the units are consistent.

"""

function load(T,fnames::Vector{<:AbstractString},datanames::Vector{<:AbstractString};
              qv_flags = [GOOD_VALUE,PROBABLY_GOOD_VALUE],
              nametype = :P01,
              qvlocalname = "QV:SEADATANET")

    profiles = T[]
    lons = T[]
    lats = T[]
    times = DateTime[]
    depths = T[]
    ids = String[]


    for fname in fnames

        sheet = readODVspreadsheet(fname);
        sheet_P01names = listSDNparams(sheet)

        # loop over all parameters
        for dataname in datanames
            if nametype == :P01
                if !(dataname in sheet_P01names)
                    # ignore this file
                    #@show sheet_P01names
                    warn("no data in $(fname)")
                    continue
                end
            elseif nametype == :localname
                if !(dataname in localnames(sheet))
                    # ignore this file
                    #@show localnames(sheet)
                    warn("no data in $(fname)")
                    continue
                end
            end

            info("Starting loop on the $(nprofiles(sheet)) profiles")
            for iprofile = 1:nprofiles(sheet)
                    data,data_qv,obslon,obslat,obsdepth,obsdepth_qv,obstime,
                       obstime_qv,EDMO,LOCAL_CDI_ID = loadprofile(T,sheet,iprofile,dataname;
                                                                  nametype = nametype,
                                                                  qvlocalname = qvlocalname)

                    # concatenate EDMO and LOCAL_CDI_ID separated by a hypthen
                    obsids = String[e * "-" * l for (e,l) in zip(EDMO,LOCAL_CDI_ID)]

                    # select data matching quality flags

                    good_data = falses(size(data_qv))
                    for flag in qv_flags
                        good_data[:] =  good_data .| (data_qv .== flag)
                    end

                    good_time = goodflag(obstime_qv,qv_flags)
                    good_depth = goodflag(obsdepth_qv,qv_flags)

                    good = good_data .& good_time .& good_depth

                    append!(profiles,data[good])
                    append!(lons,obslon[good])
                    append!(lats,obslat[good])
                    append!(depths,obsdepth[good])
                    append!(times,obstime[good])
                    append!(ids,obsids[good])
            end
            info("Done reading the profiles")
        end
    end

    return profiles,lons,lats,depths,times,ids
end

"""
     profiles,lons,lats,depths,times,ids = load(T,dir,P01names)

Load all ODV files under the directory `dir` corresponding the
one of the parameter names `P01names`. The resulting vectors have the data
type `T` (expect `times` and `ids` which are vectors of `DateTime` and
`String` respectively).

No checks are done if the units are consistent.
"""


function load(T,dir::AbstractString,datanames;
              qv_flags = [GOOD_VALUE,PROBABLY_GOOD_VALUE],
              nametype = :P01)
    fnames = cat(1,[[joinpath(root, file) for file in files if endswith(file,".txt")] for (root, dirs, files) in walkdir(dir)]...)
    return load(T,fnames,datanames; qv_flags = qv_flags, nametype = nametype)
end

export readODVspreadsheet, listSDNparams, nprofiles

end
