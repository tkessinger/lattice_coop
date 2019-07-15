#!/usr/bin/env julia

## test_game_graph.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot results of parallel_game_graph.jl.
## Look for dependence on selection strength w.

using CSV, PyPlot, Statistics

runs = CSV.read("output/test_game_graph_3.csv")

for (ui, update_type) in enumerate(sort(unique(runs[:update_type])))
    for (ki, k) in enumerate(sort(unique(runs[:k])))
        for (wi, w) in enumerate(sort(unique(runs[:w])))
            tmp_p_c_traj = []
            tmp_p_d_traj = []
            tmp_p_cc_traj = []
            tmp_p_cd_traj = []
            tmp_p_dd_traj = []
            fixed = 0
            total_trials = 0
            tmp_runs = runs[(runs[:w] .== w) .& (runs[:k] .== k), :]
            for traj in eachrow(tmp_runs[:,:])
                total_trials += 1
                p_c_traj = parse.(Float64,String.(split(traj[:p_c_trajectory], r",|\[|\]")[2:end-1]))
                p_cc_traj = parse.(Float64,String.(split(traj[:p_cc_trajectory], r",|\[|\]")[2:end-1]))
                if p_c_traj[end] != 0 && p_c_traj[end] != 1
                    push!(tmp_p_c_traj, p_c_traj)
                    push!(tmp_p_d_traj, 1.0 .- p_c_traj)
                    push!(tmp_p_cc_traj, p_cc_traj)
                    push!(tmp_p_cd_traj, p_c_traj - p_cc_traj)
                    push!(tmp_p_dd_traj, 1.0 .- 2.0*(p_c_traj - p_cc_traj))
                else
                    fixed += 1
                end
            end
            println("k = $k, w = $w, update type = $update_type, number of trials = $total_trials, fixation prob = $(1.0*fixed/total_trials)")
            fig = plt.figure()
            plt.plot(mean(tmp_p_c_traj), label = L"\rho_c")
            plt.plot(mean(tmp_p_d_traj), label = L"\rho_d")
            plt.legend(loc=2)
            plt.title("k = $k, w = $w, update type = $update_type")
            plt.tight_layout()
            display(fig)

            fig = plt.figure()
            plt.plot(mean(tmp_p_cc_traj), label = L"\rho_{cc}")
            plt.plot(mean(tmp_p_cd_traj), label = L"\rho_{cd}")
            plt.plot(mean(tmp_p_dd_traj), label = L"\rho_{dd}")
            plt.legend(loc=2)
            plt.title("k = $k, w = $w, update type = $update_type")
            plt.tight_layout()
            display(fig)
        end

    end
end
    #plt.plot(tmp_freq_traj, label="cooperator frequency")
    #plt.plot(ones(length(tmp_freq_traj)) - tmp_freq_traj, label="defector frequency")
