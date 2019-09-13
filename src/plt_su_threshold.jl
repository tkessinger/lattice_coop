#!/usr/bin/env julia

## plt_w_csv.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Duplicate Su et al. (2019) figure 4.
## Look for dependence on changing the threshold value.

using Revise
using CSV, PyPlot, Statistics, MultiEdgeGame

# output info about the fixation probabilities?
verbose = false

# load simulation output as a dataframe
runs = CSV.read("output/su_fig_threshold.csv")
#runs = CSV.read("output/su_fig4c_3.csv")
# dicts to store fixation probabilities
c_fixation = Dict{Tuple{Int64, Int64, String}, Array{Float64, 1}}()
d_fixation = Dict{Tuple{Int64, Int64, String}, Array{Float64, 1}}()

c_extinction = Dict{Tuple{Int64, Int64, String}, Array{Float64, 1}}()
d_extinction = Dict{Tuple{Int64, Int64, String}, Array{Float64, 1}}()

analytical_result = Dict{Tuple{Int64, Int64, String}, Array{Float64, 1}}()

analytical_c_fixation = Dict{Tuple{Int64, Int64, String}, Array{Float64, 1}}()
analytical_d_fixation = Dict{Tuple{Int64, Int64, String}, Array{Float64, 1}}()

k = sort(unique(runs[:k]))[1]
w = sort(unique(runs[:w]))[1]
N = sort(unique(runs[:N]))[1]
#accumulate_payoffs = sort(unique(runs[:accumulate]))[1]
#accumulate_payoffs ? (game = DoL_payoff_accumulate) : (game = DoL_payoff)

# get unique values from the runs dataframe
g1_vals = sort(unique(runs[:g1]))
g1_vals = [2,3]
bc_ratio_vals = sort(unique(runs[:BC_ratio]))
thresholds = sort(unique(runs[:threshold]))
#thresholds = [2]
fitness_functions = sort(unique(runs[:fitness_function]))
#fitness_functions = [fitness_functions[1]]

num_trials = sort(unique(runs[:num_trials]))[1]

# iterate over all update types, ks, and ws
for (ti, threshold_val) in enumerate(thresholds)
    for (gi, gi_val) in enumerate(g1_vals)
        for (fi, f_val) in enumerate(fitness_functions)
            game = get_game_function(f_val)
            key = (gi_val, threshold_val, f_val)
            g1 = gi_val
            g2 = k - gi_val
            c_fixation[key] = zeros(length(bc_ratio_vals))
            d_fixation[key] = zeros(length(bc_ratio_vals))
            c_extinction[key] = zeros(length(bc_ratio_vals))
            d_extinction[key] = zeros(length(bc_ratio_vals))
            analytical_result[key] = zeros(length(bc_ratio_vals))
            for (bci, bc_ratio) in enumerate(bc_ratio_vals)
                tmp_runs = runs[(runs[:g1] .== gi_val) .& (runs[:BC_ratio] .== bc_ratio) .& (runs[:threshold] .== threshold_val) .& (runs[:fitness_function] .== f_val), :]
                #total_trials = size(tmp_runs, 1)/2*num_trials
                total_trials = [0,0]
                for (ti, traj) in enumerate(eachrow(tmp_runs[:,:]))
                    total_trials[traj[:test_strat]] += traj[:num_trials]*traj[:runs_per_graph]
                    if traj[:test_strat] == 1
                        #total_trials[1] += traj[:num_trials]*traj[:runs_per_graph]
                        c_fixation[key][bci] += traj[:fixed]
                        c_extinction[key][bci] += traj[:extinct]
                    elseif traj[:test_strat] == 2
                        #total_trials[2] += traj[:num_trials]*traj[:runs_per_graph]
                        d_fixation[key][bci] += traj[:fixed]
                        d_extinction[key][bci] += traj[:extinct]
                    end
                end
                println("$(c_fixation[key][bci]), $(d_fixation[key][bci]), $total_trials")
                c_fixation[key][bci] /= total_trials[1]
                c_extinction[key][bci] /= total_trials[1]
                d_fixation[key][bci] /= total_trials[2]
                d_extinction[key][bci] /= total_trials[2]
                for s1 in 0:g1
                    for s2 in 0:g2
                        sigma = sigma_coeff([s1, s2], [g1, g2])
                        a = game(1, 1.0*bc_ratio, 1.0, s1, s2, threshold_val)
                        b = game(2, 1.0*bc_ratio, 1.0, (g1 - s1), (g2 - s2), threshold_val)
                        #println("$g1, $g2, $s1, $s2, $sigma, $a, $b")
                        analytical_result[key][bci] += w*sigma*(a - b)
                    end
                end
            # c_fixation[gi_val] /= (0.5*total_trials)
            # d_fixation[gi_val] /= (0.5*total_trials)
            end
        end
    end
end

color_list = ["tab:red", "tab:green", "tab:blue", "tab:pink", "tab:cyan", "tab:purple"]

for (fi, f_val) in enumerate(fitness_functions)
    fig = plt.figure(figsize = (5,5))
    for (gi, gi_val) in enumerate(g1_vals)
        for (ti, threshold_val) in enumerate(thresholds)
            key = (gi_val, threshold_val, f_val)
            c_index = gi+(ti-1)*length(g1_vals)
            plt.plot(bc_ratio_vals, c_fixation[key] - d_fixation[key],
               c = color_list[c_index], ls = "", marker = ".", label = "g1 = $gi_val, threshold = $threshold_val")
            plt.plot(bc_ratio_vals, analytical_result[key],
                c = color_list[c_index])
        end
    end
    plt.xlim([minimum(bc_ratio_vals), maximum(bc_ratio_vals)])
    plt.xlabel(L"B/C")
    plt.ylim([-w, 2*w])
    plt.title("N = $N, w = $w, fitness function = $f_val")
    plt.ylabel(L"\rho_A - \rho_B")
    plt.legend(loc=2)
    plt.tight_layout()
    display(fig)
end



# fig = plt.figure(figsize = (6,6))
# for (gi, gi_val) in enumerate(g1_vals)
#     plt.plot(bc_ratio_vals, d_fixation[gi_val],
#        c = color_list[gi], ls = "", marker = ".", label = "g1 = $gi_val, g2 = $(k - gi_val)")
#     # plt.plot(bc_ratio_vals, analytical_result[gi_val],
#     #     c = color_list[gi])
# end
# plt.xlim([minimum(bc_ratio_vals), maximum(bc_ratio_vals)])
# plt.xlabel(L"B/C")
# #plt.ylim([-0.01, 0.01])
# plt.ylabel(L"\rho_B")
# plt.legend(loc=2)
# plt.tight_layout()
# display(fig)
#
# fig = plt.figure(figsize = (6,6))
# for (gi, gi_val) in enumerate(g1_vals)
#     plt.plot(bc_ratio_vals, c_fixation[gi_val],
#        c = color_list[gi], ls = "", marker = ".", label = "g1 = $gi_val, g2 = $(k - gi_val)")
#     # plt.plot(bc_ratio_vals, analytical_result[gi_val],
#     #     c = color_list[gi])
# end
# plt.xlim([minimum(bc_ratio_vals), maximum(bc_ratio_vals)])
# plt.xlabel(L"B/C")
# #plt.ylim([-0.01, 0.01])
# plt.ylabel(L"\rho_A")
# plt.legend(loc=2)
# plt.tight_layout()
# display(fig)
