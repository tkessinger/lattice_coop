#!/usr/bin/env julia

## MultiEdgeGame.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Module for simulating multiplayer games on random regular graphs.

module MultiEdgeGame

    using Random, StatsBase

    export MultiEdgeGraph
    export evolve!, initialize_fitnesses!
    export get_frequencies, get_frequency, get_pair_frequencies, get_conditional_frequencies

    mutable struct Individual
        strategy::String
        neighbors::Array{Array{Int64, 1}}
        all_neighbors::Array{Int64, 1}
        fitness::Float64
    end

    function isregular(
        neighbor_list::Array{Array{Int64,1}},
        degree::Int64
        )
        num_appearances = zeros(Int64, length(neighbor_list))
        #println("$neighbor_list")
        for (li, list) in enumerate(neighbor_list)
            if length(list) != degree
                return false
            elseif li ∈ list
                return false
            end
            for ii in list
                num_appearances[ii] += 1
            end
        end
        #println("$num_appearances")
        if any(num_appearances .!= degree)
            return false
        else
            return true
        end
    end

    function random_graph_edges(
        N::Int64,
        k::Int64,
        gi::Array{Int64, 1}
        )

        if sum(gi) < k
            push!(gi, k - sum(gi))
        end

        all_gi_neighbors = [[Int64[] for x in 1:N] for i in 1:gi]
        all_neighbors = [Int64 for x in 1:N]

        used_edges = Array{Int64,1}[]

        for (ii, gi_neighbors) in enumerate(all_gi_neighbors)

            while !isregular(g2_neighbors, g2)
                S = collect(combinations(1:N,2))
                filter!(x -> !(x ∈ used_edges), S)
                [gi_neighbors[x] = Int64[] for x in 1:N]

                while !isempty(S)
                    #println(S)
                    wts = zeros(length(S))
                    for (pi, (i,j)) in enumerate(S)
                        wts[pi] = (gi[ii] - length(gi_neighbors[i]))*
                            (gi[ii] - length(gi_neighbors[j]))
                    end

                    #newpair = S[sample(1:length(S), Weights(wts))]
                    newpair = sample(S, Weights(wts))
                    push!(gi_neighbors[newpair[1]], newpair[2])
                    push!(gi_neighbors[newpair[2]], newpair[1])
                    filter!(x -> x != newpair, S)
                    for (pi, (i,j)) in enumerate(S)
                        if length(gi_neighbors[i]) > (gi[ii] - 1) ||
                            length(gi_neighbors[j]) > (gi[ii] - 1)
                            filter!(x -> x != [i,j], S)
                        end
                    end
                end
            end
            all_neighbors = [vcat(all_neighbors[x], gi_neighbors[x]) for x in 1:N]

            for i in 1:N
                for j in 1:gi[ii]
                    edge = [i, gi_neighbors[i][j]]
                    if issorted(edge)
                        push!(used_edges, edge)
                    end
                end
            end
            all_neighbors = g1_neighbors
        end
    end

    mutable struct MultiEdgeGraph
        # by design this is a random regular graph

        n::Int64 # number of nodes
        k::Int64 # degree of graph
        gi::Array{Int64, 1} # degrees of subgraphs, or, number of edges of each type
        w::Float64 # strength of selection
        game::Array{Float64, 2} # an associated game to be played on the graph

        # the next three are the only mutable attributes
        generation::Int64 # current generation
        pop::Array{Individual} # the population and all its members


        function GameGraph(
            n::Int64,
            k::Int64,
            w::Float64,
            game::Array{Float64, 2},
            graph::SimpleGraph
            )
            return new(n, k, w, game, graph,
                0, ones(Int64, n), zeros(n))
        end
    end

    function initialize_fitnesses!(
        gg::GameGraph
        )
        # assign the correct fitness values to each individual in the population
        # ideally this would be done within a constructor
        for indv in vertices(gg.graph)
            payoff = 0
            for neighbor in neighbors(gg.graph, indv)
                payoff += gg.game[gg.strategies[indv], gg.strategies[neighbor]]
            end
            fitness = 1 - gg.w + gg.w * payoff
            gg.fitnesses[indv] = fitness
        end
    end

    function update_indv!(
        gg::GameGraph,
        invader::Int64,
        invadee::Int64
        )
        # if an individual gets replaced, we need to update their fitness
        # as well as that of all their neighbors
        update_strategies!(gg, invader, invadee)
        update_fitness!(gg, invadee)
    end

    function update_strategies!(
        gg::GameGraph,
        invader::Int64,
        invadee::Int64
        )
        # replaces the strategy of a replaced individual
        gg.strategies[invadee] = gg.strategies[invader]
    end

    function update_fitness!(
        gg::GameGraph,
        indv::Int64
        )
        # updates the fitness of an individual whose strategy is replaced
        # as well as the fitnesses of all their neighbors
        payoff = 0.0
        fitness = 0.0
        for neighbor in neighbors(gg.graph, indv)
            payoff += gg.game[gg.strategies[indv], gg.strategies[neighbor]]
        end
        fitness = 1.0 - gg.w + gg.w * payoff
        return fitness
    end

    function get_frequencies(
        gg::GameGraph
        )
        # returns the frequencies of each strategy in the population
        freqs = zeros(size(gg.game)[1])
        for strategy in gg.strategies
            freqs[strategy] += 1
        end
        freqs /= gg.n
        return freqs
    end

    function get_frequency(
        gg::GameGraph,
        strategy::Int64
        )
        # returns the frequency of a single strategy in the population
        return count(i->(i==strategy), gg.strategies)/gg.n
    end

    function get_pair_frequencies(
        gg::GameGraph
        )
        # returns the pair frequencies
        # useful for comparison with the pair approximation
        pair_freqs = zeros(size(gg.game))
        for indv in vertices(gg.graph)
            for neighbor in neighbors(gg.graph, indv)
                pair_freqs[gg.strategies[indv], gg.strategies[neighbor]] += 1.0
            end
        end
        # pair frequencies should be symmetric: x_{ij} = x_{ji}
        pair_freqs += transpose(pair_freqs)
        # there are nk/2 total edges
        pair_freqs /= (2*gg.n*gg.k)
        return pair_freqs
    end

    function get_conditional_frequencies(
        gg::GameGraph,
        return_all::Bool=false
        )
        # returns the conditional frequencies q_{i|j} = x_{ij}/x_j
        conditional_freqs = zeros(size(gg.game))
        pair_freqs = get_pair_frequencies(gg)
        freqs = get_frequencies(gg)
        conditional_freqs = pair_freqs./transpose(freqs)
        if return_all
            return (freqs, pair_freqs, conditional_freqs)
        else
            return conditional_freqs
        end
    end

    function evolve!(
        gg::GameGraph,
        num_gens::Int64=1,
        update_rule::String="birth_death"
        )
        # depending on the update rule and their fitnesses,
        # replaces an individual with a neighbor

        for i in 1:num_gens
            if update_rule == "birth_death" || update_rule == "bd"
                # an individual is chosen to reproduce proportional to fitness
                # they replace a random neighbor
                indv = sample(1:gg.n, Weights(gg.fitnesses))
                neighbor_to_replace = rand(neighbors(gg.graph, indv))
                # this if block is added in case there's no need to update
                if gg.strategies[indv] != gg.strategies[neighbor_to_replace]
                    update_indv!(gg, indv, neighbor_to_replace)
                end
            elseif update_rule == "death_birth" || update_rule == "db"
                # an individual is chosen to die at random
                # their neighbors compete to replace them, proportional to fitness
                indv = rand(1:gg.n)
                neighbor_fitnesses = gg.fitnesses[neighbors(gg.graph, indv)]
                invading_neighbor = sample(neighbors(gg.graph, indv), Weights(neighbor_fitnesses))
                if gg.strategies[invading_neighbor] != gg.strategies[indv]
                    update_indv!(gg, invading_neighbor, indv)
                end
            elseif update_rule == "imitation" || update_rule == "im"
                # an individual compares its fitness against its neighbors and itself
                # then adopts their strategy, proportional to fitness
                indv = rand(1:gg.n)
                neighbor_fitnesses = gg.fitnesses[vcat(indv, neighbors(gg.graph, indv))]
                invading_neighbor = sample(vcat(indv, neighbors(gg.graph, indv)), Weights(neighbor_fitnesses))
                if gg.strategies[invading_neighbor] != gg.strategies[indv]
                    update_indv!(gg, invading_neighbor, indv)
                end
            end
        gg.generation += 1
        end
    end

end
