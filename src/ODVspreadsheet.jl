module ODVspreadsheet

using Logging
using StringEncodings

# Set logging level(DEBUG, INFO, WARNING, ERROR or CRITICAL)
loglevel = WARNING
Logging.configure(level=loglevel);

"""
Define composite type that will contain:
* the metadata (dictionary),
* the column labels (array) and
* the profiles (array of arrays).
"""

global Spreadsheet
type Spreadsheet
    metadata::Dict{String,String}
    columnLabels::Array{SubString{String},1}
    profileList::Array{Any,1}
end

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

    # Context manager
    open(datafile, enc"Latin1", "r") do f
        line = readline(f)

        # Read the metadata (lines starting with //)
        while line[1:2] == "//"

            # Identify metadata fields using regex
            # (name of the field is between < > and </ >)
            m = match(r"<(\w+)>(.+)</(\w+)>", line)

            if m != nothing
                debug("Match found")
                debug(m[1] * ": " * m[2])
                # Add key - value in the dictionnary
                metadata[String(m[1])] = String(m[2])
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
                for ii in nonempty_ind
                    push!(profile[ii], line[ii]);
                end
            end

        end

        # Add the last profile to the list
        push!(profileList, profile);

        info("No. of profiles in the file: " * string(nprofiles))
        ODVdata = Spreadsheet(metadata, columnLabels, profileList)
        return ODVdata
    end
end

export readODVspreadsheet

end
