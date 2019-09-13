#!/usr/bin/env julia

## plt_w_csv.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot results of parallel_game_graph.jl.
## Look for dependence on selection strength w.

using Revise
using CSV, PyPlot, Statistics, MultiEdgeGame

# output info about the fixation probabilities?
verbose = false

# load simulation output as a dataframe
#runs = CSV.read("output/su_fig_multiply_v_large_low_w.csv")
runs = CSV.read("output/su_fig4a.csv")
# dicts to store fixation probabilities
c_fixation = Dict{Int64, Array{Float64, 1}}()
d_fixation = Dict{Int64, Array{Float64, 1}}()

c_extinction = Dict{Int64, Array{Float64, 1}}()
d_extinction = Dict{Int64, Array{Float64, 1}}()

analytical_result = Dict{Int64, Array{Float64, 1}}()

analytical_c_fixation = Dict{Int64, Array{Float64, 1}}()
analytical_d_fixation = Dict{Int64, Array{Float64, 1}}()

k = sort(unique(runs[:k]))[1]
w = sort(unique(runs[:w]))[1]
N = sort(unique(runs[:N]))[1]
#accumulate_payoffs = sort(unique(runs[:accumulate]))[1]
#accumulate_payoffs ? (game = DoL_payoff_accumulate) : (game = DoL_payoff)

game = get_game_function("DoL_payoff_add")

# get unique values from the runs dataframe
g1_vals = sort(unique(runs[:g1]))
bc_ratio_vals = sort(unique(runs[:BC_ratio]))

num_trials = sort(unique(runs[:num_trials]))[1]

# iterate over all update types, ks, and ws
for (gi, gi_val) in enumerate(g1_vals)
    g1 = gi_val
    g2 = k - gi_val
    c_fixation[gi_val] = zeros(length(bc_ratio_vals))
    d_fixation[gi_val] = zeros(length(bc_ratio_vals))
    c_extinction[gi_val] = zeros(length(bc_ratio_vals))
    d_extinction[gi_val] = zeros(length(bc_ratio_vals))
    analytical_result[gi_val] = zeros(length(bc_ratio_vals))
    for (bci, bc_ratio) in enumerate(bc_ratio_vals)
        tmp_runs = runs[(runs[:g1] .== gi_val) .& (runs[:BC_ratio] .== bc_ratio), :]
        #total_trials = size(tmp_runs, 1)/2*num_trials
        total_trials = [0,0]
        for (ti, traj) in enumerate(eachrow(tmp_runs[:,:]))
            total_trials[traj[:test_strat]] += traj[:num_trials]*traj[:runs_per_graph]
            if traj[:test_strat] == 1
                #total_trials[1] += traj[:num_trials]*traj[:runs_per_graph]
                c_fixation[gi_val][bci] += traj[:fixed]
                c_extinction[gi_val][bci] += traj[:extinct]
            elseif traj[:test_strat] == 2
                #total_trials[2] += traj[:num_trials]*traj[:runs_per_graph]
                d_fixation[gi_val][bci] += traj[:fixed]
                d_extinction[gi_val][bci] += traj[:extinct]
            end
        end
        println("$(c_fixation[gi_val][bci]), $(d_fixation[gi_val][bci]), $total_trials")
        c_fixation[gi_val][bci] /= total_trials[1]
        c_extinction[gi_val][bci] /= total_trials[1]
        d_fixation[gi_val][bci] /= total_trials[2]
        d_extinction[gi_val][bci] /= total_trials[2]
        for s1 in 0:g1
            for s2 in 0:g2
                sigma = sigma_coeff([s1, s2], [g1, g2])
                a = game(1, 1.0*bc_ratio, 1.0, s1, s2)
                b = game(2, 1.0*bc_ratio, 1.0, (g1 - s1), (g2 - s2))
                #println("$g1, $g2, $s1, $s2, $sigma, $a, $b")
                analytical_result[gi_val][bci] += w*sigma*(a - b)
            end
        end
    # c_fixation[gi_val] /= (0.5*total_trials)
    # d_fixation[gi_val] /= (0.5*total_trials)
    end
end

color_list = ["tab:red", "tab:green", "tab:blue"]

fig = plt.figure(figsize = (5,5))
for (gi, gi_val) in enumerate(g1_vals)
    plt.plot(bc_ratio_vals, c_fixation[gi_val] - d_fixation[gi_val],
       c = color_list[gi], ls = "", marker = ".", label = "g1 = $gi_val, g2 = $(k - gi_val)")
    plt.plot(bc_ratio_vals, analytical_result[gi_val],
        c = color_list[gi])
end
plt.xlim([minimum(bc_ratio_vals), maximum(bc_ratio_vals)])
plt.xlabel(L"B/C")
plt.ylim([-w, 2*w])
plt.title("N = $N, w = $w")
plt.ylabel(L"\rho_A - \rho_B")
plt.legend(loc=2)
plt.tight_layout()
display(fig)

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

#[println("$keyt, $(mean(p_c_trajectories[(keyt)])[2])") for keyt in Base.product(k_vals, w_vals, update_types)]
#
# color_list = ["tab:blue", "tab:orange", "tab:green", "tab:red",
#     "tab:purple", "tab:gray", "tab:pink"]
#
# # make a 2x3 subplot
# fig, axs = plt.subplots(length(update_types), length(k_vals),
#     figsize = (14,7), sharey="row", sharex="col")
# for (ui, update_type) in enumerate(update_types)
#     for (ki, k) in enumerate(k_vals)
#         ax = axs[ui,ki]
#         for (wi, w) in enumerate(w_vals)
#             keyt = (k, w, update_type)
#             ax.plot(mean(p_c_trajectories[keyt]),
#                 label = L"w = "*string(w), c = color_list[wi])
#         end
#         # set axis limits
#         ax.set_ylim([0.3, 0.7])
#         ax.set_xlim([0, num_gens])
#         ax.set_title("k = $k, update type = $update_type")
#     end
# end
# axs[1].legend(loc=2)
# # set column labels
# axs[1,1].set_ylabel(L"\rho_c")
# axs[2,1].set_ylabel(L"\rho_c")
# # set row labels
# [axs[2,ii].set_xlabel("time") for ii in 1:length(k_vals)]
# plt.tight_layout()
# display(fig)
# plt.savefig("figures/ohtsuki_trajectories.pdf")
#
# # make a 2x3 subplot
# fig, axs = plt.subplots(length(update_types), length(k_vals),
#     figsize = (14,7), sharey="row", sharex="col")
# for (ui, update_type) in enumerate(update_types)
#     for (ki, k) in enumerate(k_vals)
#         ax = axs[ui,ki]
#         for (wi, w) in enumerate(w_vals)
#             keyt = (k, w, update_type)
#             ax.plot(mean(p_cc_trajectories[keyt]),
#                 label = L"w = "*string(w), c = color_list[wi])
#             #ax.plot(1.0 .- mean(p_cc_trajectories[keyt]) - 2*(mean(p_c_trajectories[keyt]) .- mean(p_cc_trajectories[keyt])),
#             #    label = L"\rho_{dd}, w = "*string(w), c = color_list[wi], ls = "dotted")
#         end
#         # set axis limits
#         ax.set_ylim([0.0, 1.0])
#         ax.set_xlim([0, num_gens])
#         ax.set_title("k = $k, update type = $update_type")
#     end
# end
# axs[1].legend(loc=2)
# # set column labels
# axs[1,1].set_ylabel(L"\rho_{cc}")
# axs[2,1].set_ylabel(L"\rho_{cc}")
# # set row labels
# [axs[2,ii].set_xlabel("time") for ii in 1:length(k_vals)]
# plt.tight_layout()
# display(fig)
# plt.savefig("figures/ohtsuki_p_cc_trajectories.pdf")
#
# # make a 2x3 subplot
# fig, axs = plt.subplots(length(update_types), length(k_vals),
#     figsize = (14,7), sharey="row", sharex="col")
# for (ui, update_type) in enumerate(update_types)
#     for (ki, k) in enumerate(k_vals)
#         ax = axs[ui,ki]
#         for (wi, w) in enumerate(w_vals)
#             keyt = (k, w, update_type)
#             ax.plot(1.0 .- mean(p_cc_trajectories[keyt]) .-
#                 2*(mean(p_c_trajectories[keyt]) .-
#                 mean(p_cc_trajectories[keyt])),
#                 label = L"w = "*string(w), c = color_list[wi])
#         end
#         # set axis limits
#         ax.set_ylim([0.0, 1.0])
#         ax.set_xlim([0, num_gens])
#         ax.set_title("k = $k, update type = $update_type")
#     end
# end
# axs[1].legend(loc=2)
# # set column labels
# axs[1,1].set_ylabel(L"\rho_{dd}")
# axs[2,1].set_ylabel(L"\rho_{dd}")
# # set row labels
# [axs[2,ii].set_xlabel("time") for ii in 1:length(k_vals)]
# plt.tight_layout()
# display(fig)
# plt.savefig("figures/ohtsuki_p_dd_trajectories.pdf")
#
# fig, axs = plt.subplots(length(update_types), length(k_vals),
#     figsize = (14,7), sharey="row", sharex="col")
# for (ui, update_type) in enumerate(update_types)
#     for (ki, k) in enumerate(k_vals)
#         ax = axs[ui,ki]
#         for (wi, w) in enumerate(w_vals)
#             keyt = (k, w, update_type)
#             ax.plot(mean(p_c_trajectories[keyt]) .-
#                 mean(p_cc_trajectories[keyt]),
#                 label = L"w = "*string(w), c = color_list[wi])
#         end
#         # set axis limits
#         ax.set_ylim([0.0, 1.0])
#         ax.set_xlim([0, num_gens])
#         ax.set_title("k = $k, update type = $update_type")
#     end
# end
# axs[1].legend(loc=2)
# # set column labels
# axs[1,1].set_ylabel(L"\rho_{cd}")
# axs[2,1].set_ylabel(L"\rho_{cd}")
# # set row labels
# [axs[2,ii].set_xlabel("time") for ii in 1:length(k_vals)]
# plt.tight_layout()
# display(fig)
# plt.savefig("figures/ohtsuki_p_cd_trajectories.pdf")
#
# # make a 1x2 subplot
# fig, axs = plt.subplots(1, length(update_types),
#     figsize = (8,4), sharey="all")
# for (ui, update_type) in enumerate(sort(unique(runs[:update_type])))
#     ax = axs[ui]
#     for (ki, k) in enumerate(sort(unique(runs[:k])))
#         ax.plot(w_vals, [c_fixation[(k, w, update_type)] for w in w_vals],
#             label = "c, k = $k", c = color_list[ki])
#         ax.plot(w_vals, [d_fixation[(k, w, update_type)] for w in w_vals],
#             label = "d, k = $k", c = color_list[ki], ls = "dotted")
#     end
#     # set axis limits, labels, and titles
#     ax.set_xscale("log")
#     ax.set_ylim([0, 1])
#     ax.set_xlim([minimum(w_vals), maximum(w_vals)])
#     ax.set_xlabel(L"w")
#     ax.set_title("update type = $update_type")
# end
# axs[1].legend(loc=2)
# # set row label
# axs[1].set_ylabel("fixation probability")
# plt.tight_layout()
# display(fig)
# plt.savefig("figures/ohtsuki_fixation_probs.pdf")
