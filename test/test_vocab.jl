using Base.Test
import divand



#collection = divand.Vocab.Collection("http://vocab.nerc.ac.uk/collection/P01/current/")

collection = divand.Vocab.SDNCollection("P01")


url = "http://vocab.nerc.ac.uk/collection/P01/current/PSALPR01/"
collect,tag,key = divand.Vocab.splitURL(url)

@test collect == "P01"
@test tag == "current"
@test key == "PSALPR01"

concept = collection["PSALPR01"]
@test contains(divand.Vocab.prefLabel(concept),"salinity")
@test contains(divand.Vocab.notation(concept),"P01")
@test contains(divand.Vocab.altLabel(concept),"sal")
@test typeof(divand.Vocab.date(concept)) == DateTime
