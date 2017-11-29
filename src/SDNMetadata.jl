"""encode parameters as key-value separated by : and +"""
encodeWMSStyle(params) = join([k * ':' * string(v) for (k,v) in  params ],"+")


"""
    ncglobalattrib,ncvarattrib = SDNMetadata(metadata,fi)

Based on the information in the dictionary `metadata` and the analysed 4D field
`fi` produce a list of NetCDF global and variable attributes for `divand_save2`.
"""

function SDNMetadata(metadata,filename,varname,lonr,latr,fi)

    const layersep = "*"

    pathname = joinpath(dirname(@__FILE__),"..")

    const PROJECTS = Dict(
        "EMODNET-chemistry" =>  Dict(
            "baseurl_visualization" => "http://ec.oceanbrowser.net/emodnet/",
            "baseurl_wms" => "http://ec.oceanbrowser.net/emodnet/Python/web/wms",
            "baseurl_http" => "http://ec.oceanbrowser.net/data/emodnet-domains",
            "baseurl_opendap" =>  "http://ec.oceanbrowser.net:8081/data/emodnet-domains",
            "template" => joinpath(pathname,"templates","emodnet-chemistry.xml"),
        ),
        "SeaDataNet" => Dict(
            "baseurl_visualization" => "http://sdn.oceanbrowser.net/web-vis/",
            "baseurl_wms" => "http://sdn.oceanbrowser.net/web-vis/Python/web/wms",
            "baseurl_http" => "http://sdn.oceanbrowser.net/data/SeaDataNet-domains",
            "baseurl_opendap" => "http://sdn.oceanbrowser.net:8081/data/SeaDataNet-domains",
            "template" => joinpath(pathname,"templates","seadatanet.xml"),
        ),
        "SeaDataCloud" => Dict(
            "baseurl_visualization" => "http://sdn.oceanbrowser.net/web-vis/",
            "baseurl_wms" => "http://sdn.oceanbrowser.net/web-vis/Python/web/wms",
            "baseurl_http" => "http://sdn.oceanbrowser.net/data/SeaDataNet-domains",
            "baseurl_opendap" => "http://sdn.oceanbrowser.net:8081/data/SeaDataNet-domains",
            "template" => joinpath(pathname,"templates","seadatanet.xml"),
        )
    )



    sdn_parameter_urn = metadata["parameter_keyword_urn"]
    sdn_parameter = Vocab.resolve(sdn_parameter_urn)
    sdn_parameter_name = Vocab.prefLabel(sdn_parameter)

    # units
    sdn_uom = Vocab.findfirst(sdn_parameter,"related","P06")
    sdn_uom_urn = Vocab.urn(sdn_uom)
    sdn_uom_name = Vocab.prefLabel(sdn_uom)


    # derived attributes

    project = PROJECTS[metadata["project"]]

    ncglobalattrib = OrderedDict{String,String}()

    institution = Vocab.resolve(metadata["institution_urn"])
    area = Vocab.resolve(metadata["area_keywords_urn"][1])
    parameter = Vocab.resolve(metadata["parameter_keyword_urn"])
    product_version =  metadata["product_version"]
    product_type = "ANA"

    for (k,v) in metadata
        if k in ["netcdf_standard_name","netcdf_long_name","netcdf_units"]
            # these attribute do not need to be saved as global attributes
            continue
        end

        if k == "institution_urn"
            # Name from EDMO registry from http://seadatanet.maris2.nl/edmo/
            ncglobalattrib["institution"] = Vocab.name(institution)
        elseif k == "parameter_keyword_urn"
            # Preferred label from SeaDataNet Vocabulary P35
            # example: Water body ammonium
            ncglobalattrib["parameter_keyword"] = Vocab.prefLabel(parameter)
        elseif k == "search_keywords_urn"
            # Preferred label from SeaDataNet Parameter Discovery Vocabulary P02
            # example: Ammonium and ammonia concentration parameters in water bodies
            ncglobalattrib["search_keywords"] = join(Vocab.prefLabel.(Vocab.resolve.(v)),", ")
        elseif k == "area_keywords_urn"
            # Preferred label from SeaDataNet Vocabulary C19
            # example: Black Sea
            ncglobalattrib["area_keywords"] = join(Vocab.prefLabel.(Vocab.resolve.(v)),", ")
        elseif k == "product_version"
            # External shortname
            # example: SISMER-Atlantic Sea-Water_body_silicate-1.0-ANA
            ncglobalattrib["product_code"] = join([
                Vocab.abbreviated_name(institution),
                Vocab.prefLabel(area),
                Vocab.prefLabel(parameter),
                product_version,
                product_type],"-")
        end

        if typeof(v) <: Vector
            ncglobalattrib[k] = join(v,", ")
        else
            ncglobalattrib[k] = v
        end

    end

    # Where the product can be downloaded
    ncglobalattrib["data_access"] = project["baseurl_http"]

    # Where the product can be visualised
    ncglobalattrib["WEB_visualisation"] = project["baseurl_visualization"]

    # default layer (first time instance and surface)
    field = fi[:,:,end,1]
    default_field_min,default_field_max = extrema(field[.!isnan.(field)])

    ncglobalattrib["preview"] = string(HTTP.URL(project["baseurl_wms"]; query =
                                                OrderedDict(
                                                    "service" => "WMS",
                                                    "request" => "GetMap",
                                                    "version" => "1.3.0",
                                                    "crs" => "CRS:84",
                                                    "bbox" => "$(lonr[1]),$(latr[1]),$(lonr[end]),$(latr[end])",
                                                    "decorated" => "true",
                                                    "format" => "image/png",
                                                    "layers" => Vocab.prefLabel(area) * "/" * filename * layersep * varname,
                                                    "styles" => encodeWMSStyle(Dict("vmin" => default_field_min,
                                                                                    "vmax" => default_field_max)),
                                                    #"elevation" => "-0.0",
                                                    #"time" => "03",
                                                    "transparent" => "true",
                                                    "height" => "500",
                                                    "width" => "800")))

ncvarattrib = OrderedDict("units" => metadata["netcdf_units"],
                          "standard_name" => metadata["netcdf_standard_name"],
                          "long_name" => metadata["netcdf_long_name"],
                          #  "sdn_parameter_urn" => sdn_parameter_urn,
                          #  "sdn_parameter_name" => sdn_parameter_name,
                          #  "sdn_uom_urn" => sdn_uom_urn,
                          #  "sdn_uom_name" => sdn_uom_name,
                          )
# Example  http://ec.oceanbrowser.net/emodnet/Python/web/wms?styles=vmax%3A1.97872%2Bvmin%3A1.67568&format=image%2Fpng&height=500&bbox=27.5%2C42.0%2C30.0%2C44.5&decorated=true&transparent=true&layers=Back+Sea%2Fvarname.4Danl_autumn.nc%2Avarname&crs=CRS%3A84&service=WMS&request=GetMap&width=800&version=1.3.0

return ncglobalattrib,ncvarattrib
end
