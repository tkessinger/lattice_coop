#!/usr/bin/env julia

## plt_su_DoL_analytics.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot analytical predictions for B/C ratio from Su et al. (2019).

using Revise
using CSV, PyPlot, Statistics, MultiEdgeGame

k = 15
g1_vals = collect(1:k)
#g1_vals = collect(1:5)

threshold_vals = [1,2,3]
flat_DoL_BC_ratio_vals = Dict{Int64, Array{Float64, 1}}()
add_DoL_BC_ratio_vals = Dict{Int64, Array{Float64, 1}}()
#mult_DoL_BC_ratio_vals = zeros(Float64, k-1)
for (ti, threshold_val) in enumerate(threshold_vals)
    flat_DoL_BC_ratio_vals[threshold_val] = zeros(Float64, k)
    add_DoL_BC_ratio_vals[threshold_val] = zeros(Float64, k)
    for (g1, g1_val) in enumerate(g1_vals)
        g = [g1_val, k-g1_val]
        flat_DoL_BC_ratio_vals[threshold_val][g1] += BC_ratio_DoL(k, g, DoL_flat, 4, threshold_val)
        add_DoL_BC_ratio_vals[threshold_val][g1] += BC_ratio_DoL(k, g, DoL_add, 4, threshold_val)
        #mult_DoL_BC_ratio_vals[g1] += BC_ratio_DoL(k, g, DoL_multiply, 4, threshold)
    end
end

#labels = ["threshold", "additive payoffs"]#, "multiplicative payoffs"]
color_list = ["tab:red", "tab:pink", "tab:blue", "tab:purple", "tab:orange", "tab:grey"]

labels = ["flat, threshold = 1", "flat, threshold = 2", "flat, threshold = 3",
    "additive, threshold = 1", "additive, threshold = 2", "additive, threshold = 3"]

vals_to_plot = [flat_DoL_BC_ratio_vals[1], flat_DoL_BC_ratio_vals[2], flat_DoL_BC_ratio_vals[3],
    add_DoL_BC_ratio_vals[1], add_DoL_BC_ratio_vals[2], add_DoL_BC_ratio_vals[3]]#, mult_DoL_BC_ratio_vals]
fig = plt.figure(figsize = (6,6))
for (vi, vals) in enumerate(vals_to_plot)
    plt.plot(g1_vals, vals,
       c = color_list[vi], label=labels[vi])
end
plt.xlim([0, k])
#plt.ylim([0, 1])
plt.xlabel(L"g_1")
plt.ylabel(L"B/C*")
plt.legend(loc=2)
plt.tight_layout()
display(fig)
