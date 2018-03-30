module Vocab

using Base
using EzXML
import Requests
import Base.find
import Base.findfirst
import Base.repr

const namespaces = Dict(
    "rdf" => "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    "skos" => "http://www.w3.org/2004/02/skos/core#",
    "dc" => "http://purl.org/dc/terms/",
    "dce" => "http://purl.org/dc/elements/1.1/",
    "rdfs" => "http://www.w3.org/2000/01/rdf-schema#",
    "owl" => "http://www.w3.org/2002/07/owl#" )

const CFStandardNameURL = "http://cfconventions.org/Data/cf-standard-names/current/src/cf-standard-name-table.xml"
const CFAreaTypesURL = "http://cfconventions.org/Data/area-type-table/current/src/area-type-table.xml"

type CFVocab
    xdoc :: EzXML.Document
end

type CFEntry
    node :: EzXML.Node
end


"""
    collection = Vocab.CFVocab()
    collection = Vocab.CFVocab(url = url)

Create a Dict-like object represeting the NetCDF CF Standard Name vocabulary.
If the `url` is not provided then current CF Standard Name list
$(CFStandardNameURL) is used. Individual standard names are retirved by indexing
which return an object of the type `CFEntry`:

```julia
collection = Vocab.CFVocab()
entry = collection["sea_water_temperature"]
```


"""
function CFVocab(; url = CFStandardNameURL)
    r = Requests.get(url)
    data = readstring(r)
    try
        xdoc = EzXML.parsexml(data)
    catch
        @show data[1:100]
        error("unable to parse XML file")
    end
    return CFVocab(xdoc)
end

function Base.getindex(c::CFVocab,stdname::AbstractString)
    return CFEntry([e for e in find(c.xdoc,"//entry") if e["id"] == stdname ][1])
end

"""
    bool = haskey(collection::CFVocab,stdname)

Return true if `stdname` is part of the NetCDF CF Standard Name vocabulary
`collection`.
"""

Base.haskey(c::CFVocab,stdname) = length([e for e in find(c.xdoc,"//entry") if e["id"] == stdname ]) > 0

for (method,tag) in [(:description,"description"),
                     (:canonical_units,"canonical_units")]
    @eval begin
"""
    str = description(entry::CFEntry)
    str = canonical_units(entry::CFEntry)

    Return the description or the canonical units of the `entry`.
"""
        $method(e::CFEntry) = nodecontent(findfirst(e.node,$tag))
    end
end

"""
    collection,tag,key = Vocab.splitURL(url)

Split a concept URL into collection, tag and key.
url must finishe with a slash.
"""

function splitURL(url)
    collection,tag,key = split(url,"/")[end-3:end-1]
    return collection,tag,key
end


"""
    entry = Vocab.resolve(urn)

Resolve a SeaDataNet URN (Uniform Resource Name) and returns the corresponding
EDMO entry or Vocabulary concept. For example:

```julia
concept = Vocab.resolve("SDN:P021:current:TEMP")
```
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

#            s = Vocab.$($tag)(urn::AbstractString)
for (method,tag,docname) in [(:prefLabel,"prefLabel","preferred label"),
                             (:notation,"notation","identifier"),
                             (:altLabel,"altLabel","alternative label"),
                             (:definition,"definition","definition")]
    @eval begin
        @doc """
            s = Vocab.$($tag)(c::Vocab.Concept)

        Return the $($docname) of a concept `c`
        """ $method(c::Concept) = nodecontent(findfirst(root(c.xdoc),"//skos:" * $tag,namespaces))

        @doc """
            s = Vocab.$($tag)(urn::AbstractString)

        Return the $($docname) of a concept usings it URN (Uniform Resource Name)
        """ $method(urn::AbstractString) = $method(resolve(urn))
    end
end

URL(c::Concept) = findfirst(root(c.xdoc),"//skos:Concept",namespaces)["rdf:about"]

date(c::Concept) = DateTime(
   nodecontent(findfirst(root(c.xdoc),"//dc:date",namespaces)),
   "yyyy-mm-dd HH:MM:SS.s")

urn(c::Concept) = notation(c)

"""
    find(c::Concept,name,collection)

Return a list of related concepts in the collection `collection`.
`name` can be the string "related", "narrower", "broader".
"""
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

"""
    findfirst(c::Concept,name,collection)

Return the first related concepts in the collection `collection`.
`name` can be the string "related", "narrower", "broader".
"""

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
