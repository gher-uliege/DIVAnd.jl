
divand


# divand.jl documentation

```@docs
divandrun
divandgo
divand_averaged_bg
load_mask
domain
SDNMetadata
divand.divand_save2
```

# Vocabulary

urn_str


```@docs
Vocab.CFVocab
haskey(collection::Vocab.CFVocab,stdname)
Vocab.SDNCollection
Vocab.prefLabel
Vocab.altLabel
Vocab.notation
Vocab.find(c::Vocab.Concept,name,collection)
Vocab.description
Vocab.canonical_units
Vocab.splitURL
```

# Information for developpers

## Update the documentation

Install

```bash
pip3 install --user mkdocs
pip3 install --user python-markdown-math
```

```julia
Pkg.add("Documenter")
```
