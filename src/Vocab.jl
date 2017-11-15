module Vocab

using Base
using EzXML
using Requests

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

type Concept
    xdoc :: EzXML.Document
end

function Concept(url::AbstractString)
    r = Requests.get(url)
    xdoc = parsexml(readstring(r))
    return Concept(xdoc)
end

prefLabel(c::Concept) = nodecontent(findfirst(root(c.xdoc),"//skos:prefLabel",namespaces))
notation(c::Concept) = nodecontent(findfirst(root(c.xdoc),"//skos:notation",namespaces))
altLabel(c::Concept) = nodecontent(findfirst(root(c.xdoc),"//skos:altLabel",namespaces))
date(c::Concept) = DateTime(
   nodecontent(findfirst(root(c.xdoc),"//dc:date",namespaces)),
   "yyyy-mm-dd HH:MM:SS.s")

# baseurl e.g. http://vocab.nerc.ac.uk/collection/P01/current/
# must containt a trailing /
type Collection{T <: AbstractString}
    baseurl :: T
end

Base.getindex(c::Collection,identifier::AbstractString) = Concept(c.baseurl * identifier)

SDNCollection(name) = Collection("http://vocab.nerc.ac.uk/collection/$(name)/current/")

end
