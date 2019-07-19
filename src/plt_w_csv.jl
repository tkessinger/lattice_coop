#!/usr/bin/env julia

## plt_w_csv.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot results of parallel_game_graph.jl.
## Look for dependence on selection strength w.

using CSV, PyPlot, Statistics

# output info about the fixation probabilities?
verbose = false

# load simulation output as a dataframe
runs = CSV.read("output/test_game_graph_3.csv")

# dicts to store fixation probabilities
c_fixation = Dict{Tuple{Float64, Float64, String},Float64}()
d_fixation = Dict{Tuple{Float64, Float64, String},Float64}()

# get unique values from the runs dataframe
w_vals = sort(unique(runs[:w]))
k_vals = sort(unique(runs[:k]))
update_types = sort(unique(runs[:update_type]))

# this will be needed for plotting
num_gens = maximum(runs[:num_gens])

# initialize the trajectory array dicts
function initialize_traj_dicts(
    )
    return Dict{Tuple{Float64, Float64, String},Array{Array{Float64,1}}}()
end
p_c_trajectories = initialize_traj_dicts()
# note that p_d = 1.0 - p_c
p_cc_trajectories = initialize_traj_dicts()
# note that p_cd = p_c - p_cc
# and p_dd = 1.0 .- 2.0(p_cd)

# iterate over all update types, ks, and ws
for (ui, update_type) in enumerate(update_types)
    for (ki, k) in enumerate(k_vals)
        for (wi, w) in enumerate(w_vals)
            # declare key
            keyt = (k, w, update_type)
            # initialize trajectory arrays
            # they will be arrays of arrays
            p_c_trajectories[keyt] = []
            p_cc_trajectories[keyt] = []
            c_fixed = 0
            d_fixed = 0
            total_trials = 0
            # get all trajectories satisfying w, k, and update_type
            tmp_runs = runs[(runs[:w] .== w) .& (runs[:update_type] .== update_type) .& (runs[:k] .== k), :]
            for (ti, traj) in enumerate(eachrow(tmp_runs[:,:]))
                total_trials += 1
                # trajectories are currently saved as strings
                # this reads them as arrays
                p_c_traj = parse.(Float64,String.(split(traj[:p_c_trajectory], r",|\[|\]")[2:end-1]))
                p_cc_traj = parse.(Float64,String.(split(traj[:p_cc_trajectory], r",|\[|\]")[2:end-1]))
                # if nothing fixed, keep the trajectory
                if p_c_traj[end] != 0 && p_c_traj[end] != 1
                    push!(p_c_trajectories[keyt], p_c_traj)
                    push!(p_cc_trajectories[keyt], p_cc_traj)
                    # push!(tmp_p_c_traj, p_c_traj)
                    # push!(tmp_p_d_traj, 1.0 .- p_c_traj)
                    # push!(tmp_p_cc_traj, p_cc_traj)
                    # push!(tmp_p_cd_traj, p_c_traj - p_cc_traj)
                    # push!(tmp_p_dd_traj, 1.0 .- 2.0*(p_c_traj - p_cc_traj))
                # else, figure out what fixed
                elseif p_c_traj[end] == 0
                    d_fixed += 1
                elseif p_c_traj[end] == 1
                    c_fixed += 1
                end
            end
            # compute fixation probabilities
            c_fixation[keyt] = (1.0*c_fixed/total_trials)
            d_fixation[keyt] = (1.0*d_fixed/total_trials)
            if verbose
                printstring = "k = $k, w = $w, update type = $update_type, "
                printstring *= "number of trials = $total_trials, cooperator fixation = $(1.0*c_fixed/total_trials), defector fixation = $(1.0*d_fixed/total_trials)"
                println(printstring)
            end
        end

    end
end

#[println("$keyt, $(mean(p_c_trajectories[(keyt)])[2])") for keyt in Base.product(k_vals, w_vals, update_types)]

color_list = ["tab:blue", "tab:orange", "tab:green", "tab:red",
    "tab:purple", "tab:gray", "tab:pink"]

# make a 2x3 subplot
fig, axs = plt.subplots(length(update_types), length(k_vals),
    figsize = (14,7), sharey="row", sharex="col")
for (ui, update_type) in enumerate(update_types)
    for (ki, k) in enumerate(k_vals)
        ax = axs[ui,ki]
        for (wi, w) in enumerate(w_vals)
            keyt = (k, w, update_type)
            ax.plot(mean(p_c_trajectories[keyt]),
                label = L"w = "*string(w), c = color_list[wi])
        end
        # set axis limits
        ax.set_ylim([0.3, 0.7])
        ax.set_xlim([0, num_gens])
        ax.set_title("k = $k, update type = $update_type")
    end
end
axs[1].legend(loc=2)
# set column labels
axs[1,1].set_ylabel(L"\rho_c")
axs[2,1].set_ylabel(L"\rho_c")
# set row labels
[axs[2,ii].set_xlabel("time") for ii in 1:length(k_vals)]
plt.tight_layout()
display(fig)
plt.savefig("figures/ohtsuki_trajectories.pdf")

# make a 2x3 subplot
fig, axs = plt.subplots(length(update_types), length(k_vals),
    figsize = (14,7), sharey="row", sharex="col")
for (ui, update_type) in enumerate(update_types)
    for (ki, k) in enumerate(k_vals)
        ax = axs[ui,ki]
        for (wi, w) in enumerate(w_vals)
            keyt = (k, w, update_type)
            ax.plot(mean(p_cc_trajectories[keyt]),
                label = L"w = "*string(w), c = color_list[wi])
            #ax.plot(1.0 .- mean(p_cc_trajectories[keyt]) - 2*(mean(p_c_trajectories[keyt]) .- mean(p_cc_trajectories[keyt])),
            #    label = L"\rho_{dd}, w = "*string(w), c = color_list[wi], ls = "dotted")
        end
        # set axis limits
        ax.set_ylim([0.0, 1.0])
        ax.set_xlim([0, num_gens])
        ax.set_title("k = $k, update type = $update_type")
    end
end
axs[1].legend(loc=2)
# set column labels
axs[1,1].set_ylabel(L"\rho_{cc}")
axs[2,1].set_ylabel(L"\rho_{cc}")
# set row labels
[axs[2,ii].set_xlabel("time") for ii in 1:length(k_vals)]
plt.tight_layout()
display(fig)
plt.savefig("figures/ohtsuki_p_cc_trajectories.pdf")

# make a 2x3 subplot
fig, axs = plt.subplots(length(update_types), length(k_vals),
    figsize = (14,7), sharey="row", sharex="col")
for (ui, update_type) in enumerate(update_types)
    for (ki, k) in enumerate(k_vals)
        ax = axs[ui,ki]
        for (wi, w) in enumerate(w_vals)
            keyt = (k, w, update_type)
            ax.plot(1.0 .- mean(p_cc_trajectories[keyt]) .-
                2*(mean(p_c_trajectories[keyt]) .-
                mean(p_cc_trajectories[keyt])),
                label = L"w = "*string(w), c = color_list[wi])
        end
        # set axis limits
        ax.set_ylim([0.0, 1.0])
        ax.set_xlim([0, num_gens])
        ax.set_title("k = $k, update type = $update_type")
    end
end
axs[1].legend(loc=2)
# set column labels
axs[1,1].set_ylabel(L"\rho_{dd}")
axs[2,1].set_ylabel(L"\rho_{dd}")
# set row labels
[axs[2,ii].set_xlabel("time") for ii in 1:length(k_vals)]
plt.tight_layout()
display(fig)
plt.savefig("figures/ohtsuki_p_dd_trajectories.pdf")

fig, axs = plt.subplots(length(update_types), length(k_vals),
    figsize = (14,7), sharey="row", sharex="col")
for (ui, update_type) in enumerate(update_types)
    for (ki, k) in enumerate(k_vals)
        ax = axs[ui,ki]
        for (wi, w) in enumerate(w_vals)
            keyt = (k, w, update_type)
            ax.plot(mean(p_c_trajectories[keyt]) .-
                mean(p_cc_trajectories[keyt]),
                label = L"w = "*string(w), c = color_list[wi])
        end
        # set axis limits
        ax.set_ylim([0.0, 1.0])
        ax.set_xlim([0, num_gens])
        ax.set_title("k = $k, update type = $update_type")
    end
end
axs[1].legend(loc=2)
# set column labels
axs[1,1].set_ylabel(L"\rho_{cd}")
axs[2,1].set_ylabel(L"\rho_{cd}")
# set row labels
[axs[2,ii].set_xlabel("time") for ii in 1:length(k_vals)]
plt.tight_layout()
display(fig)
plt.savefig("figures/ohtsuki_p_cd_trajectories.pdf")

# make a 1x2 subplot
fig, axs = plt.subplots(1, length(update_types),
    figsize = (8,4), sharey="all")
for (ui, update_type) in enumerate(sort(unique(runs[:update_type])))
    ax = axs[ui]
    for (ki, k) in enumerate(sort(unique(runs[:k])))
        ax.plot(w_vals, [c_fixation[(k, w, update_type)] for w in w_vals],
            label = "c, k = $k", c = color_list[ki])
        ax.plot(w_vals, [d_fixation[(k, w, update_type)] for w in w_vals],
            label = "d, k = $k", c = color_list[ki], ls = "dotted")
    end
    # set axis limits, labels, and titles
    ax.set_xscale("log")
    ax.set_ylim([0, 1])
    ax.set_xlim([minimum(w_vals), maximum(w_vals)])
    ax.set_xlabel(L"w")
    ax.set_title("update type = $update_type")
end
axs[1].legend(loc=2)
# set row label
axs[1].set_ylabel("fixation probability")
plt.tight_layout()
display(fig)
plt.savefig("figures/ohtsuki_fixation_probs.pdf")
