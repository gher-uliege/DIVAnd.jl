using Documenter
using DIVAnd

CI = get(ENV, "CI", nothing) == "true"

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

if CI
    deploydocs(
        repo = "github.com/gher-ulg/DIVAnd.jl.git",
    )
end
