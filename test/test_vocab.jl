using Base.Test
import divand


# Test CF names

collection = divand.Vocab.CFVocab()
@test haskey(collection,"sea_water_temperature")

entry = collection["sea_water_temperature"]

@test contains(divand.Vocab.description(entry),"water")
@test divand.Vocab.canonical_units(entry) == "K"


#collection = divand.Vocab.Collection("http://www.seadatanet.org/urnurl/collection/P01/current/")

collection = divand.Vocab.SDNCollection("P01")


url = "http://www.seadatanet.org/urnurl/collection/P01/current/PSALPR01/"
collectionname,tag,key = divand.Vocab.splitURL(url)

@test collectionname == "P01"
@test tag == "current"
@test key == "PSALPR01"

concept = collection["PSALPR01"]
@test contains(divand.Vocab.prefLabel(concept),"salinity")
@test contains(divand.Vocab.notation(concept),"P01")
@test contains(divand.Vocab.altLabel(concept),"sal")
@test typeof(divand.Vocab.date(concept)) == DateTime


edmo = divand.Vocab.EDMO()
entry = edmo[1495]
@test typeof(divand.Vocab.name(entry)) == String
@test typeof(divand.Vocab.phone(entry)) == String
@test typeof(divand.Vocab.address(entry))  == String
@test typeof(divand.Vocab.city(entry)) == String
@test typeof(divand.Vocab.zipcode(entry)) == String
@test typeof(divand.Vocab.email(entry)) == String
@test typeof(divand.Vocab.website(entry)) == String
@test typeof(divand.Vocab.fax(entry)) == String
@test typeof(divand.Vocab.country(entry)) == String


collection = divand.Vocab.SDNCollection("P35")
concept = collection["WATERTEMP"]

units = lowercase(divand.Vocab.prefLabel(divand.Vocab.findfirst(concept,"related","P06")))

@test units == "degrees celsius"

label = divand.Vocab.prefLabel(divand.Vocab.resolve("SDN:P021:current:TEMP"))
@test label == "Temperature of the water column"

label = divand.Vocab.prefLabel(divand.Vocab.resolve("SDN:P021::TEMP"))
@test label == "Temperature of the water column"

name = divand.Vocab.name(divand.Vocab.resolve("SDN:EDMO::575"))
@test name == "National Oceanographic Data Committee"

name = divand.Vocab.name(divand.Vocab.urn"SDN:EDMO::575")
@test name == "National Oceanographic Data Committee"


