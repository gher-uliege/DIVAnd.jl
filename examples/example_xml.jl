

using Base.Test
import divand


basedir = joinpath(dirname(@__FILE__),"..","..",
                   "divand-example-data","CDI-list")

ignore_errors = true
filepath = joinpath(basedir,"Water_body_Salinity.4Danl.nc")
varname = "Salinity"
project = "SeaDataCloud"
#project = "EMODNET-chemistry"

xmlfilename = "test.xml"

# file CDI-list-export.zip is available at:
# download("http://emodnet-chemistry.maris2.nl/download/export.zip","CDI-list-export.zip")

cdilist = joinpath(basedir,"CDI-list-export.zip")
cdilist = joinpath(basedir,"export.csv")

divand.divadoxml(filepath,varname,project,cdilist,xmlfilename,
          ignore_errors = ignore_errors)
