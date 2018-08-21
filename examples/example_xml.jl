

if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
using DIVAnd


basedir = joinpath(dirname(@__FILE__),"..","..",
                   "DIVAnd-example-data","CDI-list")

ignore_errors = true
filepath = joinpath(basedir,"Water_body_Salinity.4Danl.nc")
varname = "Salinity"
project = "SeaDataCloud"
#project = "EMODNET-chemistry"

xmlfilename = "test.xml"

cdilist = joinpath(basedir,"CDI-list-export.zip")

# using the unzipped file is usually faster
#cdilist = joinpath(basedir,"export.csv")

if !isfile(cdilist)
    # file CDI-list-export.zip is available at:
    download("http://emodnet-chemistry.maris2.nl/download/export.zip",cdilist)
end


DIVAnd.divadoxml(filepath,varname,project,cdilist,xmlfilename,
          ignore_errors = ignore_errors)
