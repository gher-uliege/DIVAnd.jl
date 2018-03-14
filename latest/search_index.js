var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "divand.jl documentation",
    "title": "divand.jl documentation",
    "category": "page",
    "text": "divand"
},

{
    "location": "index.html#divand.jl-documentation-1",
    "page": "divand.jl documentation",
    "title": "divand.jl documentation",
    "category": "section",
    "text": "divand.divandrun\ndivand.divandgo\ndivand.divand_averaged_bg\ndivand.load_mask\ndivand.domain\ndivand.SDNMetadata\ndivand.save\ndivand.saveobs\ndivand.loadobs\ndivand.loadbigfile\ndivand.checkobs\ndivand.fit_isotropic\ndivand.fit\ndivand.smoothfilter\ndivand.Anam.loglin\ndivand.Anam.logit\ndivand.divadoxml\ndivand.random\ndivand.distance\ndivand.metric\ndivand.interp\ndivand.backgroundfile"
},

{
    "location": "index.html#Examples-1",
    "page": "divand.jl documentation",
    "title": "Examples",
    "category": "section",
    "text": "To run the example, you need to install PyPlot. In the folder examples of divand, you can run e.g. the example divand_simple_example_1D.jl by issuing:# cd(\"/path/to/divand/examples\")\ninclude(\"divand_simple_example_1D.jl\")Replace /path/to/divand/ by the installation directory of divand which is the output of Pkg.dir(\"divand\") if you installed divand using Julias package manager."
},

{
    "location": "index.html#divand.Vocab.CFVocab",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.CFVocab",
    "category": "type",
    "text": "collection = Vocab.CFVocab()\ncollection = Vocab.CFVocab(url = url)\n\nCreate a Dict-like object represeting the NetCDF CF Standard Name vocabulary. If the url is not provided then current CF Standard Name list http://cfconventions.org/Data/cf-standard-names/current/src/cf-standard-name-table.xml is used. Individual standard names are retirved by indexing which return an object of the type CFEntry:\n\ncollection = Vocab.CFVocab()\nentry = collection[\"sea_water_temperature\"]\n\n\n\n"
},

{
    "location": "index.html#Base.haskey-Tuple{divand.Vocab.CFVocab,Any}",
    "page": "divand.jl documentation",
    "title": "Base.haskey",
    "category": "method",
    "text": "bool = haskey(collection::CFVocab,stdname)\n\nReturn true if stdname is part of the NetCDF CF Standard Name vocabulary collection.\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.SDNCollection",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.SDNCollection",
    "category": "function",
    "text": "collection = SDNCollection(name)\n\nOpen the SeaDataNet collection with the name name at the URL http://www.seadatanet.org/urnurl/collection/ The collection can be indexed with brackets using the identifier.\n\nusing divand\ncollection = Vocab.SDNCollection(\"P01\")\nconcept = collection[\"PSALPR01\"]\n@show Vocab.prefLabel(concept)\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.prefLabel",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.prefLabel",
    "category": "function",
    "text": "s = Vocab.prefLabel(c::Vocab.Concept)\n\nReturn the preferred label of a concept c\n\n\n\ns = Vocab.prefLabel(urn::AbstractString)\n\nReturn the preferred label of a concept usings it URN (Uniform Resource Name)\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.altLabel",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.altLabel",
    "category": "function",
    "text": "s = Vocab.altLabel(c::Vocab.Concept)\n\nReturn the alternative label of a concept c\n\n\n\ns = Vocab.altLabel(urn::AbstractString)\n\nReturn the alternative label of a concept usings it URN (Uniform Resource Name)\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.notation",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.notation",
    "category": "function",
    "text": "s = Vocab.notation(c::Vocab.Concept)\n\nReturn the identifier of a concept c\n\n\n\ns = Vocab.notation(urn::AbstractString)\n\nReturn the identifier of a concept usings it URN (Uniform Resource Name)\n\n\n\n"
},

{
    "location": "index.html#Base.find-Tuple{divand.Vocab.Concept,Any,Any}",
    "page": "divand.jl documentation",
    "title": "Base.find",
    "category": "method",
    "text": "find(c::Concept,name,collection)\n\nReturn a list of related concepts in the collection collection. name can be the string \"related\", \"narrower\", \"broader\".\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.description",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.description",
    "category": "function",
    "text": "str = description(entry::CFEntry)\nstr = canonical_units(entry::CFEntry)\n\nReturn the description or the canonical units of the `entry`.\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.canonical_units",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.canonical_units",
    "category": "function",
    "text": "str = description(entry::CFEntry)\nstr = canonical_units(entry::CFEntry)\n\nReturn the description or the canonical units of the `entry`.\n\n\n\n"
},

{
    "location": "index.html#divand.Vocab.splitURL",
    "page": "divand.jl documentation",
    "title": "divand.Vocab.splitURL",
    "category": "function",
    "text": "collection,tag,key = Vocab.splitURL(url)\n\nSplit a concept URL into collection, tag and key. url must finishe with a slash.\n\n\n\n"
},

{
    "location": "index.html#Vocabulary-1",
    "page": "divand.jl documentation",
    "title": "Vocabulary",
    "category": "section",
    "text": "urn_strVocab.CFVocab\nhaskey(collection::Vocab.CFVocab,stdname)\nVocab.SDNCollection\nVocab.prefLabel\nVocab.altLabel\nVocab.notation\nVocab.find(c::Vocab.Concept,name,collection)\nVocab.description\nVocab.canonical_units\nVocab.splitURL"
},

{
    "location": "index.html#Information-for-developers-1",
    "page": "divand.jl documentation",
    "title": "Information for developers",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Update-the-documentation-1",
    "page": "divand.jl documentation",
    "title": "Update the documentation",
    "category": "section",
    "text": "InstallPkg.add(\"Documenter\")"
},

{
    "location": "index.html#Troubleshooting-1",
    "page": "divand.jl documentation",
    "title": "Troubleshooting",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#No-plot-windows-1",
    "page": "divand.jl documentation",
    "title": "No plot windows",
    "category": "section",
    "text": "If the following command doesn\'t produce any figureusing PyPlot\nplot(0, 1)a possible solution is to modify the backend: this is done by editing the python configuration file matplotlibrc. The location of this file is obtained in python with:import matplotlib\nmatplotlib.matplotlib_fnamewhich, in my case, returns \'~/.config/matplotlib/matplotlibrc\'"
},

]}
