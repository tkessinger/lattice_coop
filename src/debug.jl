#!/usr/bin/env julia

## lattice_coop.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Track cooperator frequency as a function of transaction cost.

using Random
using PyPlot
using Revise
using LaTeXStrings

using LatticeCoop

b = 3.0 # benefit to defecting
c = 1.0 # cost for initiating transaction
κ = 0.000001 # temperature parameter

verbose = true # turn this on for detailed error tracking

N = 3
numgens = 1

histories = []
overall_freqs = []

init_freq = 0.5

pop = LatticePopulation(b, c, κ, N, init_freq, true)

history = []
push!(history, pop.lattice)
for i in 1:numgens
    evolve!(pop)
    push!(history, pop.lattice)
end
