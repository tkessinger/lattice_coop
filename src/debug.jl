#!/usr/bin/env julia

## debug.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Debug code for LatticeCoop.jl.

using Random
using PyPlot
using Revise
using LaTeXStrings

using LatticeCoop

b = 3.0 # benefit to defecting
c = 1.0 # cost for initiating transaction
κ = 0.000001 # temperature parameter

verbose = true # track errors

N = 2
numgens = 1

histories = []
overall_freqs = []

init_freq = 0.5

pop = LatticePopulation(b, c, κ, N, init_freq, verbose)

history = []
push!(history, pop.lattice)
for i in 1:numgens
    evolve!(pop)
    push!(history, pop.lattice)
end
