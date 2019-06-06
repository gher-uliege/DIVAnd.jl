using Documenter
using DIVAnd

makedocs(
    format = Documenter.HTML(),
    modules = [DIVAnd],
    sitename = "DIVAnd",
    pages = [
        "index.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.

deploydocs(
    repo = "github.com/gher-ulg/DIVAnd.jl.git",
)
