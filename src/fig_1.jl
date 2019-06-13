#!/usr/bin/env julia

## fig_1.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Duplicates figure 1 of Li et al. (2019):
## cooperator frequency as a function of b and c.

using Random
using Statistics
using DataFrames
using LaTeXStrings
using LatticeCoop

using PyPlot

numgens = 1500
N = 100

clr = ["r", "b", "y", "g", "m"]

b_step = 0.01
c_step = 0.05

b_vals = collect(1.0:b_step:1.06)
c_vals = collect(0.0:c_step:1.0)
κ = 0.1

ρ_c_vals = zeros(length(b_vals), length(c_vals))
to_average = 500

for (bi, b) in enumerate(b_vals)
    for (ci, c) in enumerate(c_vals)
        println("evolving with benefit $b and transaction cost $c")
        pop = LatticePopulation(b, c, κ, N)
        history = []
        push!(history, pop.lattice)
        for i in 1:numgens
            evolve!(pop)
            push!(history, pop.lattice)
        end
        coop_freqs = 1.0/N^2*[sum(x) for x in history]
        ρ_c_vals[bi, ci] = mean(coop_freqs[end-to_average:end])
    end
end

fig = plt.figure()
color_list = ["red", "blue", "yellow", "green", "purple", "gray", "pink"]
for (bi, b) in enumerate(b_vals)
    plt.plot(c_vals, ρ_c_vals[bi,:], label="b = $b", c=color_list[bi])
end
ax=plt.gca()
ax.set_xlabel("c")
ax.set_ylabel(L"\rho_c")
ax.set_ylim([0.0,0.7])
ax.set_xlim([0.0,1.0])
plt.legend(loc=1)
plt.tight_layout()
display(fig)
