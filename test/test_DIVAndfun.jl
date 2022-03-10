using DIVAnd
using Statistics
using DelimitedFiles
using LinearAlgebra
using StableRNGs
using Test


rng = StableRNG(123)
myfunny=DIVAndfun((rand(rng,10),rand(rng,10)),rand(rng,10))
@test myfunny(0.5,0.5) ≈ 0.3904677899955127
myfunny=DIVAndfun((rand(rng,10),),rand(rng,10))
@test myfunny(0.5) ≈ 0.5072531965422283
myfunny=DIVAndfun((rand(rng,10),),rand(rng,10);xi=(collect(0.:0.1:1.0),))
@test myfunny(0.1) ≈ 0.24017172590061936