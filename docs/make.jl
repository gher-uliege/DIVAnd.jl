using Documenter
using divand

makedocs(
    format = :html,
    modules = [divand],
    sitename = "divand",
    pages = [
        "index.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.

deploydocs(
    repo = "github.com/gher-ulg/divand.jl.git",
    target = "build",
    julia  = "0.6",
    deps = nothing,
    make = nothing,
)
