#!/usr/bin/env julia

## MultiEdgeGame.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Module for simulating multiplayer games on random regular graphs.

module MultiEdgeGame

    using Random, StatsBase
    export MultiGameGraph, Population
    export DoL_payoff
    export initialize_fitnesses!
    export evolve!
    export get_frequencies, get_frequency
    export get_pair_frequencies, get_conditional_frequencies

    struct MultiGameGraph{GameFunc <: Function}

        n::Int64
        g::Array{Int64, 1}
        num_strategies::Int64
        neighbors::Array{Array{Array{Int64, 1}, 1}, 1}
        all_neighbors::Array{Array{Int64, 1}, 1}
        game_func::GameFunc
        w::Float64
        game_params::Array{Float64, 1}

    end

    mutable struct Population

        n::Int64
        strategies::Array{Int64, 1}
        fitnesses::Array{Float64, 1}
        generation::Int64

    end

    function initialize_fitnesses!(
        pop::Population,
        mgg::MultiGameGraph
        )
        # assign the correct fitness values to each individual in the population
        # ideally this would be done within a constructor
        for indv in 1:mgg.n
            payoff = mgg.game_func(pop, mgg, indv)
            fitness = 1 - mgg.w + mgg.w * payoff
            pop.fitnesses[indv] = fitness
        end
    end

    function DoL_payoff(
        pop::Population,
        mgg::MultiGameGraph,
        indv::Int64
        )
        B = mgg.game_params[1]
        C = mgg.game_params[2]
        s1 = sum([pop.strategies[x] .== 1 for x in mgg.neighbors[1][indv]]) # number of cooperators along edgetype 1
        s2 = sum([pop.strategies[x] .== 1 for x in mgg.neighbors[2][indv]]) # number of cooperators along edgetype 2
        # if indv is a cooperator
        if pop.strategies[indv] == 1
            if s2 >= 1
                return (s1 + s2 + 1)*(B-C)
            else
                return -C
            end
        # if indv is a defector
        elseif pop.strategies[indv] == 2
            if s1 >= 1 && s2 >= 1
                return (s1 + s2)*B
            else
                return 0
            end
        end
    end

    function get_frequencies(
          pop::Population,
          mgg::MultiGameGraph
          )
          # returns the frequencies of each strategy in the population
          freqs = zeros(mgg.num_strategies)
          for strategy in pop.strategies
              freqs[strategy] += 1
          end
          freqs /= pop.n
          return freqs
      end

      function get_frequency(
          pop::Population,
          mgg::MultiGameGraph,
          strategy::Int64
          )
          # returns the frequency of a single strategy in the population
          return count(i->(i==strategy), pop.strategies)/pop.n
      end

      function get_pair_frequencies(
          pop::Population,
          mgg::MultiGameGraph
          )
          # returns the pair frequencies
          # useful for comparison with the pair approximation
          pair_freqs = zeros(mgg.num_strategies, mgg.num_strategies)
          for indv in 1:pop.n
              for neighbor in mgg.all_neighbors[indv]
                  pair_freqs[pop.strategies[indv], pop.strategies[neighbor]] += 1.0
              end
          end
          # pair frequencies should be symmetric: x_{ij} = x_{ji}
          pair_freqs += transpose(pair_freqs)
          # there are nk/2 total edges
          pair_freqs /= (2*mgg.n*sum(mgg.g))
          return pair_freqs
      end

      function get_conditional_frequencies(
          pop::Population,
          mgg::MultiGameGraph,
          return_all::Bool=false
          )
          # returns the conditional frequencies q_{i|j} = x_{ij}/x_j
          conditional_freqs = zeros(mgg.num_strategies, mgg.num_strategies)
          pair_freqs = get_pair_frequencies(pop, mgg)
          freqs = get_frequencies(pop, mgg)
          conditional_freqs = pair_freqs./transpose(freqs)
          if return_all
              return (freqs, pair_freqs, conditional_freqs)
          else
              return conditional_freqs
          end
      end

    function update_indv!(
        pop::Population,
        mgg::MultiGameGraph,
        invader::Int64,
        invadee::Int64
        )
        # if an individual gets replaced, we need to update their fitness
        # as well as that of all their neighbors
        update_strategies!(pop, mgg, invader, invadee)
        update_fitness!(pop, mgg, invadee)
    end

    function update_strategies!(
        pop::Population,
        mgg::MultiGameGraph,
        invader::Int64,
        invadee::Int64
        )
        # replaces the strategy of a replaced individual
        pop.strategies[invadee] = pop.strategies[invader]
    end

    function update_fitness!(
        pop::Population,
        mgg::MultiGameGraph,
        indv::Int64
        )
        # updates the fitness of an individual whose strategy is replaced
        # as well as the fitnesses of all their neighbors
        payoff = mgg.game_func(pop, mgg, indv)
        fitness = 1.0 - mgg.w + mgg.w * payoff
        return fitness
    end

    function evolve!(
        pop::Population,
        mgg::MultiGameGraph,
        num_gens::Int64=1,
        update_rule::String="death_birth"
        )
        # depending on the update rule and their fitnesses,
        # replaces an individual with a neighbor

        for i in 1:num_gens
            if update_rule == "birth_death" || update_rule == "bd"
                # an individual is chosen to reproduce proportional to fitness
                # they replace a random neighbor
                indv = sample(1:pop.n, Weights(pop.fitnesses))
                neighbor_to_replace = rand(mgg.all_neighbors[indv])
                # this if block is added in case there's no need to update
                if pop.strategies[indv] != pop.strategies[neighbor_to_replace]
                    update_indv!(pop, mgg, indv, neighbor_to_replace)
                end
            elseif update_rule == "death_birth" || update_rule == "db"
                # an individual is chosen to die at random
                # their neighbors compete to replace them, proportional to fitness
                indv = rand(1:pop.n)
                neighbor_fitnesses = pop.fitnesses[mgg.all_neighbors[indv]]
                invading_neighbor = sample(mgg.all_neighbors[indv], Weights(neighbor_fitnesses))
                if pop.strategies[invading_neighbor] != pop.strategies[indv]
                    update_indv!(pop, mgg, invading_neighbor, indv)
                end
            elseif update_rule == "imitation" || update_rule == "im"
                # an individual compares its fitness against its neighbors and itself
                # then adopts their strategy, proportional to fitness
                indv = rand(1:pop.n)
                neighbor_fitnesses = pop.fitnesses[vcat(indv, mgg.all_neighbors[indv])]
                invading_neighbor = sample(vcat(indv, mgg.all_neighbors[indv]), Weights(neighbor_fitnesses))
                if pop.strategies[invading_neighbor] != pop.strategies[indv]
                    update_indv!(pop, mgg, invading_neighbor, indv)
                end
            end
        pop.generation += 1
        end
    end

end
#
# module MultiEdgeGame
#
#     using Random, StatsBase
#
#     export MultiGameGraph
#     export generate_multitype_connected_graph
#     export DoL_payoff, evolve!, initialize_fitnesses!
#     export get_frequencies, get_frequency, get_pair_frequencies, get_conditional_frequencies
#
#     mutable struct MultiGameGraph{Game<:Function}
#         # by design this is a random regular graph
#
#         n::Int64 # number of nodes
#         g::Array{Int64, 1} # degrees of edge-typed subgraphs
#         w::Float64 # strength of selection
#         game::Game # an associated game to be played on the graph
#         #game_params::Array{Float64, 1}
#         graph::Array{Array{Array{Int64, 1}, 1}, 1} # a Dict of individuals and their neighbors by type
#
#         # the next three are the only mutable attributes
#         generation::Int64 # current generation
#         strategies::Array{Int64, 1} # an array of individuals' strategies
#         fitnesses::Array{Float64, 1} # individuals' fitnesses depending on,
#             # e.g., neighbor payoffs
#     end
#
#     function DoL_payoff(
#         graph::Array{Array{Array{Int64, 1}, 1}, 1},
#         strategies::Array{Int64, 1},
#         indv::Int64,
#         )
#         B = 1.1
#         C = 1.1
#         s1 = sum([strategies[x] .== 1 for x in graph[indv][1]]) # number of cooperators along edgetype 1
#         s2 = sum([strategies[x] .== 1 for x in graph[indv][2]]) # number of cooperators along edgetype 2
#         # if indv is a cooperator
#         if strategies[indv] == 1
#             if s2 >= 1
#                 return (s1 + s2 + 1)*(B-C)
#             else
#                 return -C
#             end
#         # if indv is a defector
#         elseif strategies[indv] == 2
#             if s1 >= 1 && s2 >= 1
#                 return (s1 + s2)*B
#             else
#                 return 0
#             end
#         end
#     end
#
    # function initialize_fitnesses!(
    #     gg::MultiGameGraph
    #     )
    #     # assign the correct fitness values to each individual in the population
    #     # ideally this would be done within a constructor
    #     for indv in keys(gg.graph)
    #         payoff = gg.game(gg.graph, gg.strategies, indv)
    #         fitness = 1 - gg.w + gg.w * payoff
    #         gg.fitnesses[indv] = fitness
    #     end
    # end
#
#     function update_indv!(
#         gg::MultiGameGraph,
#         invader::Int64,
#         invadee::Int64
#         )
#         # if an individual gets replaced, we need to update their fitness
#         # as well as that of all their neighbors
#         update_strategies!(gg, invader, invadee)
#         update_fitness!(gg, invadee)
#     end
#
#     function update_strategies!(
#         gg::MultiGameGraph,
#         invader::Int64,
#         invadee::Int64
#         )
#         # replaces the strategy of a replaced individual
#         gg.strategies[invadee] = gg.strategies[invader]
#     end
#
#     function update_fitness!(
#         gg::MultiGameGraph,
#         indv::Int64
#         )
#         # updates the fitness of an individual whose strategy is replaced
#         # as well as the fitnesses of all their neighbors
#         payoff = gg.game(gg.graph, gg.strategies, indv)
#         fitness = 1.0 - gg.w + gg.w * payoff
#         return fitness
#     end
#
#     function get_frequencies(
#         gg::MultiGameGraph
#         )
#         # returns the frequencies of each strategy in the population
#         freqs = zeros(size(gg.game)[1])
#         for strategy in gg.strategies
#             freqs[strategy] += 1
#         end
#         freqs /= gg.n
#         return freqs
#     end
#
#     function get_frequency(
#         gg::MultiGameGraph,
#         strategy::Int64
#         )
#         # returns the frequency of a single strategy in the population
#         return count(i->(i==strategy), gg.strategies)/gg.n
#     end
#
#     function get_pair_frequencies(
#         gg::MultiGameGraph
#         )
#         # returns the pair frequencies
#         # useful for comparison with the pair approximation
#         pair_freqs = zeros(size(gg.game))
#         for indv in vertices(gg.graph)
#             for neighbor in neighbors(gg.graph, indv)
#                 pair_freqs[gg.strategies[indv], gg.strategies[neighbor]] += 1.0
#             end
#         end
#         # pair frequencies should be symmetric: x_{ij} = x_{ji}
#         pair_freqs += transpose(pair_freqs)
#         # there are nk/2 total edges
#         pair_freqs /= (2*gg.n*gg.k)
#         return pair_freqs
#     end
#
#     function get_conditional_frequencies(
#         gg::MultiGameGraph,
#         return_all::Bool=false
#         )
#         # returns the conditional frequencies q_{i|j} = x_{ij}/x_j
#         conditional_freqs = zeros(size(gg.game))
#         pair_freqs = get_pair_frequencies(gg)
#         freqs = get_frequencies(gg)
#         conditional_freqs = pair_freqs./transpose(freqs)
#         if return_all
#             return (freqs, pair_freqs, conditional_freqs)
#         else
#             return conditional_freqs
#         end
#     end
#
#     function evolve!(
#         gg::MultiGameGraph,
#         num_gens::Int64=1,
#         update_rule::String="death_birth"
#         )
#         # depending on the update rule and their fitnesses,
#         # replaces an individual with a neighbor
#
#         for i in 1:num_gens
#             if update_rule == "birth_death" || update_rule == "bd"
#                 # an individual is chosen to reproduce proportional to fitness
#                 # they replace a random neighbor
#                 indv = sample(1:gg.n, Weights(gg.fitnesses))
#                 neighbor_to_replace = rand(neighbors(gg.graph, indv))
#                 # this if block is added in case there's no need to update
#                 if gg.strategies[indv] != gg.strategies[neighbor_to_replace]
#                     update_indv!(gg, indv, neighbor_to_replace)
#                 end
#             elseif update_rule == "death_birth" || update_rule == "db"
#                 # an individual is chosen to die at random
#                 # their neighbors compete to replace them, proportional to fitness
#                 indv = rand(1:gg.n)
#                 neighbor_fitnesses = gg.fitnesses[neighbors(gg.graph, indv)]
#                 invading_neighbor = sample(neighbors(gg.graph, indv), Weights(neighbor_fitnesses))
#                 if gg.strategies[invading_neighbor] != gg.strategies[indv]
#                     update_indv!(gg, invading_neighbor, indv)
#                 end
#             elseif update_rule == "imitation" || update_rule == "im"
#                 # an individual compares its fitness against its neighbors and itself
#                 # then adopts their strategy, proportional to fitness
#                 indv = rand(1:gg.n)
#                 neighbor_fitnesses = gg.fitnesses[vcat(indv, neighbors(gg.graph, indv))]
#                 invading_neighbor = sample(vcat(indv, neighbors(gg.graph, indv)), Weights(neighbor_fitnesses))
#                 if gg.strategies[invading_neighbor] != gg.strategies[indv]
#                     update_indv!(gg, invading_neighbor, indv)
#                 end
#             end
#         gg.generation += 1
#         end
#     end
#
#
# end
