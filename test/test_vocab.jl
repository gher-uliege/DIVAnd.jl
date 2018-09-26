if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
import DIVAnd
using Compat

# Test CF names

collection = DIVAnd.Vocab.CFVocab()
@test haskey(collection,"sea_water_temperature")

entry = collection["sea_water_temperature"]

@test occursin("water",DIVAnd.Vocab.description(entry))
@test DIVAnd.Vocab.canonical_units(entry) == "K"


#collection = DIVAnd.Vocab.Collection("http://www.seadatanet.org/urnurl/collection/P01/current/")

collection = DIVAnd.Vocab.SDNCollection("P01")


url = "http://www.seadatanet.org/urnurl/collection/P01/current/PSALPR01/"
collectionname,tag,key = DIVAnd.Vocab.splitURL(url)

@test collectionname == "P01"
@test tag == "current"
@test key == "PSALPR01"

concept = collection["PSALPR01"]
@test occursin("salinity",DIVAnd.Vocab.prefLabel(concept) )
@test occursin("P01"     ,DIVAnd.Vocab.notation(concept)  )
@test occursin("sal"     ,DIVAnd.Vocab.altLabel(concept)  )
@test occursin("is"      ,DIVAnd.Vocab.definition(concept))
@test occursin(key       ,DIVAnd.Vocab.URL(concept)       )

@test typeof(DIVAnd.Vocab.date(concept)) == DateTime


# search by label

name = "P02"
search_keyword = "Chlorophyll pigment concentrations in water bodies"

collection = DIVAnd.Vocab.SDNCollection(name)
concepts = DIVAnd.Vocab.findbylabel(collection,[search_keyword])

label = DIVAnd.Vocab.prefLabel(concepts[1])
URL = DIVAnd.Vocab.URL(concepts[1])
@test label == search_keyword


edmo = DIVAnd.Vocab.EDMO()
entry = edmo[1495]
@test typeof(DIVAnd.Vocab.name(entry)) == String
@test typeof(DIVAnd.Vocab.phone(entry)) == String
@test typeof(DIVAnd.Vocab.address(entry))  == String
@test typeof(DIVAnd.Vocab.city(entry)) == String
@test typeof(DIVAnd.Vocab.zipcode(entry)) == String
@test typeof(DIVAnd.Vocab.email(entry)) == String
@test typeof(DIVAnd.Vocab.website(entry)) == String
@test typeof(DIVAnd.Vocab.fax(entry)) == String
@test typeof(DIVAnd.Vocab.country(entry)) == String


collection = DIVAnd.Vocab.SDNCollection("P35")
concept = collection["WATERTEMP"]

units = lowercase(DIVAnd.Vocab.prefLabel(DIVAnd.Vocab.findfirst(concept,"related","P06")))

@test units == "degrees celsius"

label = DIVAnd.Vocab.prefLabel(DIVAnd.Vocab.resolve("SDN:P02:current:TEMP"))
@test label == "Temperature of the water column"

label = DIVAnd.Vocab.prefLabel(DIVAnd.Vocab.resolve("SDN:P02::TEMP"))
@test label == "Temperature of the water column"

edmoname = DIVAnd.Vocab.name(DIVAnd.Vocab.resolve("SDN:EDMO::575"))
@test edmoname == "National Oceanographic Data Committee"

edmoname = DIVAnd.Vocab.name(DIVAnd.Vocab.urn"SDN:EDMO::575")
@test edmoname == "National Oceanographic Data Committee"

