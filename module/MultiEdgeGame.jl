#!/usr/bin/env julia

## MultiEdgeGame.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Module for simulating multiplayer games on random regular graphs.

module MultiEdgeGame

    using Random, StatsBase
    export MultiGameGraph, Population, GameParams
    export DoL_payoff, DoL_payoff_accumulate
    export no_DoL_payoff, no_DoL_payoff_accumulate
    export initialize_fitnesses!
    export evolve!
    export get_frequencies, get_frequency
    export get_pair_frequencies, get_conditional_frequencies
    export get_game_function, sigma_coeff
    export DoL_multiply, DoL_add, DoL_flat, BC_ratio_DoL
    export DoL_payoff_add, DoL_payoff_multiply, DoL_payoff_flat

    struct GameParams
        B::Float64
        C::Float64
        threshold::Int64
    end

    struct MultiGameGraph{GameFunc <: Function}

        n::Int64
        g::Array{Int64, 1}
        num_strategies::Int64
        neighbors::Array{Array{Array{Int64, 1}, 1}, 1}
        all_neighbors::Array{Array{Int64, 1}, 1}
        game_func::GameFunc
        w::Float64
        game_params::GameParams

    end

    mutable struct Population

        n::Int64
        strategies::Array{Int64, 1}
        fitnesses::Array{Float64, 1}
        generation::Int64

    end

    function get_game_function(funcname::String)

        game_func = Dict{String, Function}()
        #game_func["DoL_payoff_accumulate"] = DoL_payoff_accumulate
        game_func["DoL_payoff_flat"] = DoL_payoff_flat
        #game_func["no_DoL_payoff_accumulate"] = no_DoL_payoff_accumulate
        #game_func["no_DoL_payoff"] = no_DoL_payoff
        game_func["DoL_payoff_multiply"] = DoL_payoff_multiply
        game_func["DoL_payoff_add"] = DoL_payoff_add


        return game_func[funcname]
    end


    function initialize_fitnesses!(
        pop::Population,
        mgg::MultiGameGraph
        )
        # assign the correct fitness values to each individual in the population
        # ideally this would be done within a constructor
        for indv in 1:mgg.n
            payoff = obtain_payoff(pop, mgg, indv)
            fitness = 1 - mgg.w + mgg.w * payoff
            pop.fitnesses[indv] = fitness
        end
    end

    function obtain_payoff(
        pop::Population,
        mgg::MultiGameGraph,
        indv::Int64
        )
        strategy = pop.strategies[indv]
        B = mgg.game_params.B
        C = mgg.game_params.C
        threshold = mgg.game_params.threshold
        # number of cooperators along edgetype 1
        s1 = sum([pop.strategies[x] .== 1 for x in mgg.neighbors[1][indv]])
        # number of cooperators along edgetype 2
        s2 = sum([pop.strategies[x] .== 1 for x in mgg.neighbors[2][indv]])
        return mgg.game_func(strategy, B, C, s1, s2, threshold)
    end

    function DoL_multiply(
        s::Array{Int64, 1}
        )
        return prod(s)
    end

    function DoL_add(
        s::Array{Int64, 1}
        )
        return sum(s)
    end

    function DoL_flat(
        s::Array{Int64, 1}
        )
        return 1
    end

    function DoL_payoff_add(
        strategy::Int64,
        B::Float64,
        C::Float64,
        s1::Int64,
        s2::Int64,
        threshold::Int64=1
        )
        # division of labor payoff function with accumulative benefits
        return DoL_payoff(strategy, B, C, [s1, s2], DoL_add, threshold)
    end

    function DoL_payoff_multiply(
        strategy::Int64,
        B::Float64,
        C::Float64,
        s1::Int64,
        s2::Int64,
        threshold::Int64=1
        )
        # division of labor payoff function with accumulative benefits
        return DoL_payoff(strategy, B, C, [s1, s2], DoL_multiply, threshold)
    end

    function DoL_payoff_flat(
        strategy::Int64,
        B::Float64,
        C::Float64,
        s1::Int64,
        s2::Int64,
        threshold::Int64=1
        )
        # division of labor payoff function with accumulative benefits
        return DoL_payoff(strategy, B, C, [s1, s2], DoL_flat, threshold)
    end

    function DoL_payoff(
        strategy::Int64,
        B::Float64,
        C::Float64,
        s::Array{Int64, 1},
        accumulate_func::Function,
        threshold::Int64=1
        )
        # if indv is a cooperator
        s1, s2 = s[1], s[2]
        if strategy == 1
            if s1 >= threshold-1 && s2 >= threshold
                # if accumulate_payoffs is true, then accumulate benefits
                # otherwise, don't
                # either way, subtract the cost
                return accumulate_func([s[1]+1, s[2]])*B - C
            else
                return -C
            end
        # if indv is a defector
        elseif strategy == 2
            if s1 >= threshold && s2 >= threshold
                # if accumulate_payoffs is true, then accumulate benefits
                # otherwise, don't
                return accumulate_func(s)*B
            else
                return 0
            end
        end
    end
    #
    #
    # function DoL_payoff_accumulate(
    #     strategy::Int64,
    #     B::Float64,
    #     C::Float64,
    #     s1::Int64,
    #     s2::Int64,
    #     threshold::Int64=1
    #     )
    #     # division of labor payoff function with accumulative benefits
    #     return DoL_payoff(strategy, B, C, s1, s2, true)
    # end
    #
    # function DoL_payoff_multiply(
    #     strategy::Int64,
    #     B::Float64,
    #     C::Float64,
    #     s1::Int64,
    #     s2::Int64,
    #     threshold::Int64=1
    #     )
    #     # division of labor payoff function
    #     # multiplicative payoffs
    #     # it makes no sense to have an "accumulate_payoffs" parameter
    #     accumulate_payoffs = true
    #
    #     # if indv is a cooperator
    #     if strategy == 1
    #         if s1 >= threshold-1 && s2 >= threshold
    #             # if accumulate_payoffs is true, then accumulate benefits
    #             # otherwise, don't
    #             # either way, subtract the cost
    #             return (accumulate_payoffs ? (s1 + 1)*s2 : 1)*B - C
    #         else
    #             return -C
    #         end
    #     # if indv is a defector
    #     elseif strategy == 2
    #         if s1 >= threshold && s2 >= threshold
    #             # if accumulate_payoffs is true, then accumulate benefits
    #             # otherwise, don't
    #             return (accumulate_payoffs ? s1*s2 : 1)*B
    #         else
    #             return 0
    #         end
    #     end
    # end
    #
    # function DoL_payoff(
    #     strategy::Int64,
    #     B::Float64,
    #     C::Float64,
    #     s1::Int64,
    #     s2::Int64,
    #     accumulate_payoffs::Bool=false,
    #     threshold::Int64=1
    #     )
    #     # division of labor payoff function
    #
    #     # if indv is a cooperator
    #     if strategy == 1
    #         if s1 >= threshold-1 && s2 >= threshold
    #             # if accumulate_payoffs is true, then accumulate benefits
    #             # otherwise, don't
    #             # either way, subtract the cost
    #             return (accumulate_payoffs ? (s1 + s2 + 1) : 1)*B - C
    #         else
    #             return -C
    #         end
    #     # if indv is a defector
    #     elseif strategy == 2
    #         if s1 >= threshold && s2 >= threshold
    #             # if accumulate_payoffs is true, then accumulate benefits
    #             # otherwise, don't
    #             return (accumulate_payoffs ? (s1 + s2) : 1)*B
    #         else
    #             return 0
    #         end
    #     end
    # end

    function no_DoL_payoff_accumulate(
        strategy::Int64,
        B::Float64,
        C::Float64,
        s1::Int64,
        s2::Int64
        )
        # division of labor payoff function with accumulative benefits
        return no_DoL_payoff(strategy, B, C, s1, s2, true)
    end

    function no_DoL_payoff(
        strategy::Int64,
        B::Float64,
        C::Float64,
        s1::Int64,
        s2::Int64,
        accumulate_payoffs::Bool=false
        )
        # division of labor payoff function

        s = s1 + s2

        # if indv is a cooperator
        if strategy == 1
            if s >= 1
                # if accumulate_payoffs is true, then accumulate benefits
                # otherwise, don't
                # either way, subtract the cost
                return (accumulate_payoffs ? (s + 1) : 1)*B - C
            else
                return -C
            end
        # if indv is a defector
        elseif strategy == 2
            if s >= 2
                # if accumulate_payoffs is true, then accumulate benefits
                # otherwise, don't
                return (accumulate_payoffs ? s : 1)*B
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
          # normalize the frequencies
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
          # there are nk/2 total edges (k = sum(g))
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
          # return_all allows us to return all the relevant frequencies
          # might as well, since we had to calculate them anyway
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
        for neighb in mgg.all_neighbors[invadee]
            update_fitness!(pop, mgg, neighb)
        end
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
        # note that this _does not_ update the neighbor fitnesses
        # that has to be done manually: see update_indv! above
        payoff = obtain_payoff(pop, mgg, indv)
        pop.fitnesses[indv] = 1.0 - mgg.w + mgg.w * payoff
    end

    function evolve!(
        pop::Population,
        mgg::MultiGameGraph,
        num_gens::Int64=1,
        update_rule::String="death_birth",
        verbose::Bool=false
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
                #println("pop looks like $(pop.fitnesses) and $(pop.strategies)")
                indv = rand(1:pop.n)
                neighbor_fitnesses = [pop.fitnesses[x] for x in mgg.all_neighbors[indv]]
                #neighbor_fitnesses = pop.fitnesses[mgg.all_neighbors[indv]]
                invading_neighbor = sample(mgg.all_neighbors[indv], Weights(neighbor_fitnesses))

                if verbose
                    println("pop looks like $(pop.fitnesses) and $(pop.strategies)")
                    println("$indv was chosen to die")
                    println("$indv's neighbors are $(mgg.neighbors[1][indv]) and $(mgg.neighbors[2][indv])")
                    println("$indv's neighbors look like $(pop.strategies[mgg.neighbors[1][indv]]), $(pop.strategies[mgg.neighbors[2][indv]])")
                    println("neighbor_fitnesses is $neighbor_fitnesses")
                    println("$invading_neighbor was chosen to invade")
                end
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

    function Psi(
        k::Int64,
        i::Int64,
        l::Int64
        )
        return binomial(l, k-1-i)/((k-2)*(k-1)^l) + binomial(k-1-l, k-i)/(1.0*k-1)^(k-1-l)
    end

    function Phi(
        k::Int64,
        i::Int64,
        l::Int64
        )
        return binomial(l, k-i)/(k-1)^l + binomial(k-1-l, k-1-i)/((k-2)*(1.0*k-1)^(k-1-l))
    end

    function sigma_coeff(
        s::Array{Int64, 1},
        g::Array{Int64, 1}
        )
        k = sum(g)
        if length(s) != length(g)
            println("error: s and g must be the same length")
        end
        first_term = (k-2)^(k - sum(s))/(k^2*(k+1)*(k+2))
        second_term = prod([binomial(g[j], s[j]) for j in 1:length(s)])/binomial(k, sum(s))
        third_term = 0
        for l in 0:k
            third_term += (k-l)*((2*k + (k-2)*l)*Psi(k, sum(s), l) + (k^2 - (k-2)*l)*Phi(k, sum(s), l))
        end
        return first_term * second_term * third_term
    end

    function BC_ratio_DoL(
        k::Int64,
        g::Array{Int64, 1},
        payoff_func::Function,
        which_term::Int64,
        threshold::Int64 = 1
        )
        first_term = 0
        second_term = 0
        for s1 in (threshold-1):g[1]
            for s2 in threshold:g[2]
                payoff = payoff_func([s1 + 1, s2])
                first_term += sigma_coeff([s1, s2], g)*payoff
                println("$s1, $s2, $(sigma_coeff([s1, s2], g)), $payoff")
            end
        end
        for s1 in 0:(g[1]-threshold)
            for s2 in 0:(g[2]-threshold)
                second_term += sigma_coeff([s1, s2], g)*payoff_func([g[1]-s1, g[2]-s2])
            end
        end
        if which_term == 1
            return first_term
        elseif which_term == 2
            return second_term
        elseif which_term == 3
            return first_term - second_term
        else
            return 1.0/(first_term - second_term)
        end
    end

end
