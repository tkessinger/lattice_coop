#!/usr/bin/env julia

## test_game_graph.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Implement and test GraphGame.

using Random
using Revise
using MultiEdgeGame, GenGraph
#using PyPlot
using Profile, ProfileView
using Dates


#rc("font", size="large")

function main(
    test_strat::Int64=1,
    verbose::Bool=false
    )

    n = 100
    g = [3,0]
    w = 0.01

    update_type = "db"

    accumulate_payoffs = false

    accumulate_payoffs ? (game = DoL_payoff_accumulate) : (game = DoL_payoff)

    num_runs = 1

    B = 0.25
    C = 1.0

    #fixed = zeros(2)
    #extinct = zeros(2)

    #for i in 1:num_runs
    graph = generate_multitype_connected_graph(n, g)
    if verbose
        println("graph looks like $(graph[1])")
    end
    #game = DoL_payoff
    pop = Population(n, ones(Int64, n)*(2 - (test_strat == 2)), zeros(n), 0)
    mgg = MultiGameGraph(n, g, 2, graph[1], graph[2], DoL_payoff, w, [B,C])
    indvs_to_flip = randperm(n)[1:n÷2]
    pop.strategies[indvs_to_flip] = ones(Int64, n÷2)*test_strat
    #pop.strategies[rand(1:n)] = test_strat
    initialize_fitnesses!(pop, mgg)

    freq = get_frequency(pop, mgg, test_strat)

    while freq ∉ (0,1)
        evolve!(pop, mgg, n, update_type, verbose)
        freq = get_frequency(pop, mgg, test_strat)
    end
    if freq == 1
        return 1
    else
        return 0
    end
#    return (fixed, extinct)
end

function main1000()
    for j in 1:2
        fixed = 0
        for i in 1:1000
            fixed += main(j, false)
        end
    println("$j, $fixed")
    end
end

#Profile.clear()
main1000()
#t1 = now()
#@profile main1000()
#t2 = now()
#Profile.print(maxdepth=8)
#println("$(t2 - t1)")

# freq_trajectory = Array{Float64, 1}[]
# pair_freq_trajectory = Array{Float64, 1}[]
# cond_freq_trajectory = Array{Float64, 1}[]



# for i in 1:num_gens
#     #println(i)
#     freq_stats = get_conditional_frequencies(pop, mgg, true)
#     push!(freq_trajectory, freq_stats[1])
#     push!(pair_freq_trajectory, vec(freq_stats[2]))
#     push!(cond_freq_trajectory, vec(freq_stats[3]))
#     evolve!(pop, mgg, n, update_type)
#     #println("$(gg.generation), $(get_frequencies(gg)), $(get_pair_frequencies(gg))")
# end
# freq_trajectory = hcat(freq_trajectory...)
# pair_freq_trajectory = hcat(pair_freq_trajectory...)
# cond_freq_trajectory = hcat(cond_freq_trajectory...)
#
# rule_labels = Dict{String, String}(
#     "bd" => "birth-death",
#     "birth_death" => "birth-death",
#     "db" => "death-birth",
#     "death_birth" => "death-birth",
#     "im" => "imitation",
#     "imitation" => "imitation"
#     )
#
# color_list = ["tab:blue", "tab:orange", "tab:green", "tab:red",
#     "tab:purple", "tab:gray", "tab:pink"]
#
#
# fig, axs = plt.subplots(1,3, figsize = (15,5), sharey="all")
#
# tmp_ax = axs[1]
# tmp_ax.plot(freq_trajectory[1,:], label="cooperator", c=color_list[1])
# tmp_ax.plot(freq_trajectory[2,:], label="defector",  c=color_list[4])
# tmp_ax.set_xlabel("time")
# tmp_ax.set_ylabel("frequency")
# tmp_ax.set_ylim([0,1])
# tmp_ax.legend(loc=2)
#
# tmp_ax = axs[2]
# tmp_ax.plot(pair_freq_trajectory[1,:], label=L"\rho_{cc}", c=color_list[1])
# tmp_ax.plot(pair_freq_trajectory[2,:], label=L"\rho_{cd}", c=color_list[2])
# tmp_ax.plot(pair_freq_trajectory[4,:], label=L"\rho_{dd}", c=color_list[4])
# tmp_ax.set_xlabel("time")
# #tmp_ax.set_ylabel("frequency")
# #tmp_ax.set_ylim([0,1])
# tmp_ax.legend(loc=2)
#
# tmp_ax = axs[3]
# tmp_ax.plot(cond_freq_trajectory[1,:], label=L"q_{c|c}", c=color_list[1])
# tmp_ax.plot(cond_freq_trajectory[2,:], label=L"q_{d|c}", c=color_list[2])
# tmp_ax.plot(cond_freq_trajectory[3,:], label=L"q_{c|d}", c=color_list[3])
# tmp_ax.plot(cond_freq_trajectory[4,:], label=L"q_{d|d}", c=color_list[4])
# tmp_ax.set_xlabel("time")
# #tmp_ax.set_ylabel("frequency")
# #tmp_ax.set_ylim([0,1])
# tmp_ax.legend(loc=2)
# fig.suptitle("w = $w, g1 = $(g[1]), g2 = $(g[2]), n = $n, rule = $(rule_labels[update_type])")
# fig.tight_layout()
# fig.subplots_adjust(top=0.93)
# display(fig)
