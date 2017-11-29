module ODVspreadsheet

using Logging
using StringEncodings

# Set logging level(DEBUG, INFO, WARNING, ERROR or CRITICAL)
loglevel = WARNING
Logging.configure(level=loglevel);

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
    profileList::Array{Any,1}
end

"""
    n = nprofiles(ODVData)

Return the number of profiles in a ODV Spreadsheet `ODVData` loaded by
`readODVspreadsheet`.
"""

nprofiles(ODVData) = length(ODVData.profileList)

function initProfileList(line)
    """
    Create an empty list of lists,
    the number of internal lists is the number of columns
    found in the ODV spreadsheet.

    Input:

    * `line`: Array{SubString{String},1} as obtained by applying `split`
            to a text line read from an ODV spreadsheet.

    Output:

    * `profile`: an Array of Arrays (one per column of ODV spreadsheet).

    List of lists is preferred because the length of each list is
    not always the same:
    * the columns storing ocean variables (e.g., temperature, depth) will
    contain several values for a given profile;
    * the columns storing information about the profile (e.g., coordinates,
    station) only have one value.
    """

    debug("Creating new profile (list of lists)")

    # Compute number of columns
    ncolumns = length(line);
    debug("No. of columns: " * string(ncolumns))

    profile = []
    for i in 1:ncolumns
        push!(profile, [line[i]])
    end

    return profile
end

function getNonEmptyInd(line)
    nonempty(x) = length(x) > 0;
    nonempty_ind = find(nonempty, line);
    return nonempty_ind;
end


function readODVspreadsheet(datafile)
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

    # metadata will be stored in a dictionary
    # ODV doc: Comment lines start with two slashes  // as first two characters
    metadata = Dict{String, String}()
    SDN_parameter_mapping = Dict{String,Dict{String,String}}()

    # Context manager
    open(datafile, enc"Latin1", "r") do f
        line = readline(f)

        # Byte Order Mark (BOM) as Latin-1
        if startswith(line,"ï»¿")
            # ignore BOM
            line = line[7:end]
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
        ncols = length(columnLabels);
        debug("No. of columns: " * string(ncols))

        # Create an array that will store all the profiles
        profileList = []

        # Loop on the lines
        nlines = 0
        profile = [];
        nprofiles = 0;

        # Read the first data line to initiate the loop
        line = split(chomp(readline(f)), "\t");
        nprofiles += 1;
        debug("Working with a header line")
        debug("Create a new, empty profile")
        profile = initProfileList(line)

        while !eof(f)
            nlines += 1;
            line = split(chomp(readline(f)), "\t");

            # Count empty values
            nonempty_ind = getNonEmptyInd(line);
            debug("Indices of the non-empty columns :")
            debug(nonempty_ind);

            # some files have only white space on the last line
            if length(nonempty_ind) == 0
                continue
            end

            # If the first value (Station) is not empty,
            # then it's a header line
            if (nonempty_ind[1] == 1)
                debug("Working with a header line")
                debug("Adding the profile to the array")
                push!(profileList, profile)

                # Initiate a profile (list of lists)
                nprofiles += 1;
                debug("Create a new, empty profile")
                profile = initProfileList(line)
            else
                debug("Adding values to the existing profile")

                # section 2.3
                # If there is no data value then the data value column is left blank with the flag field set to ‘9’

                #for ii in nonempty_ind

                # keep all
                for ii in 1:length(line)
                    push!(profile[ii], line[ii]);
                end
            end

        end

        # Add the last profile to the list
        push!(profileList, profile);

        info("No. of profiles in the file: " * string(nprofiles))
        ODVdata = Spreadsheet(metadata, SDN_parameter_mapping, columnLabels, profileList)
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
https://web.archive.org/web/20171129142108/https://www.hermetic.ch/cal_stud/chron_jdate.htm
"""

parsejd(t) = DateTime(2007,2,10) + Dates.Millisecond(round(Int,(t - 2454141.5) * (24*60*60*1000)))


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
                data[i] == fillvalue
            end
        else
            if eltype(data) <: AbstractString
                data[i] = col[i]
            else
                data[i] = parse(eltype(data),col[i])
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
    lenprof = maximum(length.(profile))

    cn_data = colnumber(sheet,locname)
    data = Vector{T}(lenprof)

    if cn_data == 0
        data[:] = fillvalue
    else
        SDNparse!(profile[cn_data],fillmode,fillvalue,data)
    end

    return data
end

"""
    data,data_qv = loaddataqv(sheet,profile,locname,fillvalue; fillmode = :repeat)

The same as `loaddata`, but now the quality flag are also loaded.
"""

function loaddataqv(sheet,profile,locname,fillvalue::T; fillmode = :repeat) where T
    locnames = localnames(sheet)
    lenprof = maximum(length.(profile))

    cn_data = colnumber(sheet,locname)
    data = Vector{T}(lenprof)
    data_qv = fill("",(lenprof,))


    if cn_data == 0
        data[:] = fillvalue
        data_qv[:] = MISSING_VALUE
    else
        SDNparse!(profile[cn_data],fillmode,fillvalue,data)

        if (cn_data < length(locnames)) && (locnames[cn_data+1] == "QV:SEADATANET")
            data_qv[:] = profile[cn_data+1]
        end
    end

    return data,data_qv
end

"""
     data,data_qv,obslon,obslat,obsdepth,obstime,obstime_qv,EDMO,LOCAL_CDI_ID =
     loadSDNprofile(T,sheet,iprofile,P01name)

Load a `iprofile`-th profile from the ODV spreadsheet `sheet` of the
parameter `P01name`. The resulting vectors have the data type `T`
(expect the quality flag and `obstime`) .
"""

function loadSDNprofile(T,sheet,iprofile,P01name)
    const fillvalue = T(NaN)
    const filldate_jd = 0.
    const filldate = parsejd(filldate_jd)

    profile = sheet.profileList[iprofile]

    locnames = localnames(sheet)
    P01names = listSDNparams(sheet)

    data_locname = localnames(sheet,P01name)[1]
    data,data_qv = loaddataqv(sheet,profile,data_locname,fillvalue)
    sz = size(data)

    #cruise = loaddata(sheet,profile,"Cruise","")
    #station = loaddata(sheet,profile,"Station","")
    #ptype = loaddata(sheet,profile,"Type","")
    EDMO = loaddata(sheet,profile,"EDMO_code","")
    LOCAL_CDI_ID = loaddata(sheet,profile,"LOCAL_CDI_ID","")

    lon = loaddata(sheet,profile,"Longitude",fillvalue)
    lat = loaddata(sheet,profile,"Latitude",fillvalue)

    depth = fill(fillvalue,sz)
    if "SDN:P01::ADEPZZ01" in P01names
        localname_depth = localnames(sheet,"SDN:P01::ADEPZZ01")
        depth[:] = loaddata(sheet,profile,localname_depth,fillvalue)
    end

    time = Vector{DateTime}(sz)
    time_qv = Vector{String}(sz)

    # chronological julian day
    if "SDN:P01::CJDY1101" in P01names
        locname_time = localnames(sheet,"SDN:P01::CJDY1101")[1]
        timedata,time_qv = loaddataqv(sheet,profile,locname_time,filldate_jd)
        time = parsejd.(timedata)
    elseif "SDN:P01::DTUT8601" in P01names
        # ISO8601 format, e.g. yyyy-mm-ddThh:mm:ss.sss

        locname_time = localnames(sheet,"SDN:P01::DTUT8601")
        time,time_qv = loaddataqv(sheet,profile,locname_time,filldate)
    else
        # hopefully not necessary
        for header in ["yyyy-mm-ddThh:mm:ss.sss",
                       "yyyy-mm-ddThh:mm:ss",
                       "yyyy-mm-ddThh:mm",
                       "yyyy-mm-ddThh",
                       "yyyy-mm-dd"]
            if header in locnames
                time,time_qv = loaddataqv(sheet,profile,header,filldate)
                break
            end
        end
    end

    return data,data_qv,lon,lat,depth,time,time_qv,EDMO,LOCAL_CDI_ID
end



"""
     profiles,lons,lats,depths,times,ids = loadSDN(T,fnames,P01names)

Load all profiles in all file from the array `fnames` corresponding the
one of the parameter names `P01names`. The resulting vectors have the data
type `T` (expect `times` and `ids` which are vectors of `DateTime` and
`String` respectively).

No checks are done if the units are consistent.
"""

function loadSDN(T,fnames,P01names; qv_flags = [GOOD_VALUE,PROBABLY_GOOD_VALUE])
    profiles = T[]
    lons = T[]
    lats = T[]
    times = DateTime[]
    depths = T[]
    ids = String[]


    for fname in fnames
        debug("Loading $(fname)")

        sheet = readODVspreadsheet(fname);
        sheet_P01names = listSDNparams(sheet)

        # loop over all parameters
        for P01name in P01names
            if P01name in sheet_P01names
                for iprofile = 1:nprofiles(sheet)

                    data,data_qv,obslon,obslat,obsdepth,obstime,
                    obstime_qv,EDMO,LOCAL_CDI_ID = loadSDNprofile(T,sheet,iprofile,P01name)

                    # concatenate EDMO and LOCAL_CDI_ID separated by a hypthen
                    obsids = String[e * "-" * l for (e,l) in zip(EDMO,LOCAL_CDI_ID)]

                    # select data matching quality flags

                    good_data = falses(size(data_qv))
                    for flag in qv_flags
                        good_data[:] =  good_data .| (data_qv .== flag)
                    end

                    good_time = falses(size(obstime_qv))
                    for flag in qv_flags
                        good_time[:] =  good_time .| (obstime_qv .== flag)
                    end


                    good = good_data .& good_time

                    append!(profiles,data[good])
                    append!(lons,obslon[good])
                    append!(lats,obslat[good])
                    append!(depths,obsdepth[good])
                    append!(times,obstime[good])
                    append!(ids,obsids[good])
                end
            end
        end
    end

    return profiles,lons,lats,depths,times,ids
end



export readODVspreadsheet, listSDNparams, nprofiles

end
