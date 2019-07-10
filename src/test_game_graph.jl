#!/usr/bin/env julia

## test_game_graph.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Implement and test GraphGame.

using Random, LightGraphs
using Revise
using GraphGame
using PyPlot



#rc("font", size="large")

R = 15
S = 0
T = 16
P = 8

n = 5000
k = 3
w = 0.01

game = Float64[R S ; T P]
graph = random_regular_graph(n, k)

update_type = "im"

gg = GameGraph(n, k, w, game, graph)
indvs_to_flip = randperm(n)[1:n÷2]
gg.strategies[indvs_to_flip] = ones(Int64, n÷2)*2
initialize_fitnesses!(gg)

num_gens = 1000

freq_trajectory = Array{Float64, 1}[]
pair_freq_trajectory = Array{Float64, 1}[]
cond_freq_trajectory = Array{Float64, 1}[]

for i in 1:num_gens
    #println(i)
    freq_stats = get_conditional_frequencies(gg, true)
    push!(freq_trajectory, freq_stats[1])
    push!(pair_freq_trajectory, vec(freq_stats[2]))
    push!(cond_freq_trajectory, vec(freq_stats[3]))
    evolve!(gg, n, update_type)
    #println("$(gg.generation), $(get_frequencies(gg)), $(get_pair_frequencies(gg))")
end
freq_trajectory = hcat(freq_trajectory...)
pair_freq_trajectory = hcat(pair_freq_trajectory...)
cond_freq_trajectory = hcat(cond_freq_trajectory...)

rule_labels = Dict{String, String}(
    "bd" => "birth-death",
    "birth_death" => "birth-death",
    "db" => "death-birth",
    "death_birth" => "death-birth",
    "im" => "imitation",
    "imitation" => "imitation"
    )

color_list = ["tab:blue", "tab:orange", "tab:green", "tab:red",
    "tab:purple", "tab:gray", "tab:pink"]


fig, axs = plt.subplots(1,3, figsize = (15,5), sharey="all")

tmp_ax = axs[1]
tmp_ax.plot(freq_trajectory[1,:], label="cooperator", c=color_list[1])
tmp_ax.plot(freq_trajectory[2,:], label="defector",  c=color_list[4])
tmp_ax.set_xlabel("time")
tmp_ax.set_ylabel("frequency")
tmp_ax.set_ylim([0,1])
tmp_ax.legend(loc=2)

tmp_ax = axs[2]
tmp_ax.plot(pair_freq_trajectory[1,:], label=L"\rho_{cc}", c=color_list[1])
tmp_ax.plot(pair_freq_trajectory[2,:], label=L"\rho_{cd}", c=color_list[2])
tmp_ax.plot(pair_freq_trajectory[4,:], label=L"\rho_{dd}", c=color_list[4])
tmp_ax.set_xlabel("time")
#tmp_ax.set_ylabel("frequency")
#tmp_ax.set_ylim([0,1])
tmp_ax.legend(loc=2)

tmp_ax = axs[3]
tmp_ax.plot(cond_freq_trajectory[1,:], label=L"q_{c|c}", c=color_list[1])
tmp_ax.plot(cond_freq_trajectory[2,:], label=L"q_{d|c}", c=color_list[2])
tmp_ax.plot(cond_freq_trajectory[3,:], label=L"q_{c|d}", c=color_list[3])
tmp_ax.plot(cond_freq_trajectory[4,:], label=L"q_{d|d}", c=color_list[4])
tmp_ax.set_xlabel("time")
#tmp_ax.set_ylabel("frequency")
#tmp_ax.set_ylim([0,1])
tmp_ax.legend(loc=2)
fig.suptitle("w = $w, k = $k, n = $n, rule = $(rule_labels[update_type])")
fig.tight_layout()
fig.subplots_adjust(top=0.93)
display(fig)

#
# plt.tight_layout()
# [ax.set(xlabel=L"c") for ax in axs[4,:]]
# [ax.set(ylabel=labels[k]) for ax in axs[:,1]]
# axs[4,1].legend(loc=2)
# fig.suptitle("i = $(strategy_string(i_strat, true)), j = $(strategy_string(j_strat, true))")
# plt.tight_layout()
# fig.subplots_adjust(top=0.93)
# display(fig)
# plt.savefig("figures/H_stats_i_" * strategy_string(i_strat) * "_j_" *
#     strategy_string(j_strat) * "_" * save_strings[k] * ".pdf")
# plt.figure()
# [plt.plot(x) for x in pair_trajectory]
#
# plt.figure()
# [plt.plot(x) for x in cond_freq_trajectory]
