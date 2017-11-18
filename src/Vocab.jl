module Vocab

using Base
using EzXML
import Requests
import HTTP
import Base.find
import Base.findfirst
import Base.repr

const namespaces = Dict(
    "rdf"=> "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    "skos" => "http://www.w3.org/2004/02/skos/core#",
    "dc" => "http://purl.org/dc/terms/",
    "dce" => "http://purl.org/dc/elements/1.1/",
    "rdfs" => "http://www.w3.org/2000/01/rdf-schema#",
    "owl" => "http://www.w3.org/2002/07/owl#" )


"""
    collection,tag,key = splitURL(url)

Split a concept URL into collection, tag and key.
url must finishe with a slash.
"""

function splitURL(url)
    collection,tag,key = split(url,"/")[end-3:end-1]
    return collection,tag,key
end


"""
SDN:P021:current:TEMP
"""
function resolve(urn)
    parts = split(urn,':')
    if parts[1] == "SDN"
        if parts[2] == "EDMO"
            return EDMOEntry("http://www.seadatanet.org/urnurl/$(urn)")
        else
            return Concept("http://www.seadatanet.org/urnurl/$(urn)")
        end
    else
       return nothing
    end
end

type Concept
    xdoc :: EzXML.Document
end

function Concept(url::AbstractString)
    r = Requests.get(url)
    xdoc = parsexml(readstring(r))
    return Concept(xdoc)
end

# """
#     s = prefLabel(c::Concept)
# Return the preferred label of a concept `c`
# """

for (method,tag) in [(:prefLabel,"prefLabel"),
                     (:notation,"notation"),
                     (:altLabel,"altLabel")]
    @eval begin
        $method(c::Concept) = nodecontent(findfirst(root(c.xdoc),"//skos:" * $tag,namespaces))
        $method(urn) = $method(resolve(urn))
    end
end

date(c::Concept) = DateTime(
   nodecontent(findfirst(root(c.xdoc),"//dc:date",namespaces)),
   "yyyy-mm-dd HH:MM:SS.s")

urn(c::Concept) = notation(c)

"""name can be related, narrower, broader"""
function Base.find(c::Concept,name,collection)
    concepts = Concept[]

    for node in find(root(c.xdoc),"//skos:" * name,namespaces)
        url = node["rdf:resource"]
        coll,tag,key = splitURL(url)
        if coll == collection
            push!(concepts,Concept(url))
        end
    end

    return concepts
end

Base.findfirst(c::Concept,name,collection) = find(c,name,collection)[1]

# baseurl e.g. http://www.seadatanet.org/urnurl/collection/P01/current/
# must containt a trailing /
type Collection{T <: AbstractString}
    baseurl :: T
end

Base.getindex(c::Collection,identifier::AbstractString) = Concept(c.baseurl * identifier)

"""
    collection = SDNCollection(name)

Open the SeaDataNet collection with the name `name` at the URL
http://www.seadatanet.org/urnurl/collection/
The collection can be indexed with brackets using the identifier.

```
using divand
collection = Vocab.SDNCollection("P01")
concept = collection["PSALPR01"]
@show Vocab.prefLabel(concept)
```
"""

SDNCollection(name) = Collection("http://www.seadatanet.org/urnurl/collection/$(name)/current/")

type EDMO{T <: AbstractString}
    baseurl :: T
end

EDMO() = EDMO("http://www.seadatanet.org/urnurl/SDN:EDMO::")

Base.getindex(e::EDMO,identifier) = EDMOEntry("$(e.baseurl)$(identifier)")


type EDMOEntry
    xdoc :: EzXML.Document
end

function EDMOEntry(url::AbstractString)
    r = Requests.get(url)
    xdoc = parsexml(readstring(r))
    return EDMOEntry(xdoc)
end


get(ee::EDMOEntry,tag) = nodecontent(findfirst(root(ee.xdoc),"//" * tag))

code(ee::EDMOEntry) = get(ee,"n_code")
name(ee::EDMOEntry) = get(ee,"name")
abbreviated_name(ee::EDMOEntry) = get(ee,"abbreviated_name")
phone(ee::EDMOEntry) = get(ee,"phone")
address(ee::EDMOEntry) = get(ee,"address")
city(ee::EDMOEntry) = get(ee,"city")
zipcode(ee::EDMOEntry) = get(ee,"zipcode")
email(ee::EDMOEntry) = get(ee,"email")
website(ee::EDMOEntry) = get(ee,"website")
fax(ee::EDMOEntry) = get(ee,"fax")
country(ee::EDMOEntry) = get(ee,"c_country")

urn(ee::EDMOEntry) = "SDN:EDMO::" * get(ee,"n_code")

repr(ee::Union{Concept,EDMOEntry}) = urn(ee)

"""
    urn"SDN:x:y:z'

Resolve a SeaDataNet URN (Uniform Resource Name) using
https://www.seadatanet.org/urnurl/
"""
macro urn_str(urn)
    resolve(urn)
end

export urn_str
end
