using Documenter
using DIVAnd

CI = get(ENV, "CI", nothing) == "true"

makedocs(
    format = Documenter.HTML(; size_threshold=1_000_000),
    modules = [DIVAnd],
    sitename = "DIVAnd",
    warnonly = true,
    pages = [
        "index.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.

if CI
    deploydocs(
        repo = "github.com/gher-uliege/DIVAnd.jl.git",
    )
end
