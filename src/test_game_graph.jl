#!/usr/bin/env julia

## test_game_graph.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Implement and test GraphGame.

using Random, LightGraphs
using Revise
using GraphGame

R = 5
S = 0
T = 8
P = 1

n = 5000
k = 3
w = 0.1

game = Float64[R S ; T P]
graph = random_regular_graph(n, k)

gg = GameGraph(n, k, w, game, graph)
indvs_to_flip = randperm(n)[1:n÷2]
gg.strategies[indvs_to_flip] = ones(Int64, n÷2)*2
initialize_fitnesses!(gg)

for i in 1:100000
    evolve!(gg)
    if (i % 100 == 0)
        println("$(get_pair_frequencies(gg))")
    end
end
