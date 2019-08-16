#!/usr/bin/env julia

## plt_su_DoL_analytics.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot analytical predictions for B/C ratio from Su et al. (2019).

using Revise
using CSV, PyPlot, Statistics, MultiEdgeGame

k = 16
g1_vals = collect(1:k-1)

flat_DoL_BC_ratio_vals = zeros(Float64, k-1)
add_DoL_BC_ratio_vals = zeros(Float64, k-1)
mult_DoL_BC_ratio_vals = zeros(Float64, k-1)
for (g1, g1_val) in enumerate(g1_vals)
    g = [g1_val, k-g1_val]
    flat_DoL_BC_ratio_vals[g1] += BC_ratio_DoL(k, g, DoL_flat, 4)
    add_DoL_BC_ratio_vals[g1] += BC_ratio_DoL(k, g, DoL_add, 4)
    mult_DoL_BC_ratio_vals[g1] += BC_ratio_DoL(k, g, DoL_multiply, 4)
end

labels = ["threshold", "additive payoffs", "multiplicative payoffs"]
color_list = ["tab:red", "tab:green", "tab:blue"]

vals_to_plot = [flat_DoL_BC_ratio_vals, add_DoL_BC_ratio_vals, mult_DoL_BC_ratio_vals]
fig = plt.figure(figsize = (6,6))
for (vi, vals) in enumerate(vals_to_plot)
    plt.plot(g1_vals, vals,
       c = color_list[vi], label=labels[vi])
end
plt.xlim([0, k])
plt.xlabel(L"g_1")
plt.ylabel(L"B/C")
plt.legend(loc=2)
plt.tight_layout()
display(fig)
