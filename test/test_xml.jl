if VERSION >= v"0.7.0-beta.0"
    using Test
else
    using Base.Test
end
import DIVAnd


db = Dict{Tuple{Int64,String},Tuple{Bool,Vector{Int64}}}(
    (1,"A") => (true,[1,2]),
    (2,"A") => (true,[10,20]),
    (2,"B") => (true,[100]),
    (2,"C") => (false,[100]),
    (2,"C") => (false,[1051]),
)

obsids = ["1-A","2-B","1000-A","2-C"]


originators,notfound =
    @static if VERSION >= v"0.7.0-beta.0"
        @test_logs (:info,r".*EDMO.*") match_mode=:any DIVAnd.get_originators_from_obsid(
            db,obsids; ignore_errors = true)
    else
        @test_warn r".*EDMO.*" DIVAnd.get_originators_from_obsid(
            db,obsids; ignore_errors = true)
    end

@test any(originator -> originator["EDMO_CODE"] == "1",originators)
@test notfound[1]["edmo"] == 1000
@test notfound[1]["local_cdi"] == "A"

contact = DIVAnd.getedmoinfo(1579,"role")
@test contact["country"] == "Belgium"


csvfile = IOBuffer(""""active","author_edmo","cdi_identifier","originator_edmo"
"True","42","A_53","12"
"True","42","A_55","13,17"
"True","42","A_57","14"
"True","42","A_59","12"
"True","42","A_61","12"
"True","42","A_63","12"
""")

db = DIVAnd.loadoriginators(csvfile)

@test db[(42,"A_53")] == (true,[12])
@test db[(42,"A_55")] == (true,[13,17])
