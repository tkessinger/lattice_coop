#!/usr/bin/env julia

## cost.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Track cooperator frequency as a function of transaction cost.

using Random
using Revise
using LaTeXStrings
using LatticeCoop

using PyPlot

b = 1.02 # benefit to defecting
κ = 0.1 # temperature parameter

histories = []
overall_freqs = []

fig = plt.figure()

numgens = 1000
N = 100

c = 0.0
for i in 1:5

    println("evolving with transaction cost $c")

    pop = LatticePopulation(b, c, κ, N)

    history = []

    push!(history, pop.lattice)
    for i in 1:numgens
        evolve!(pop)
        push!(history, pop.lattice)
    end

    push!(histories, history)

    coop_freqs = 1.0/N^2*[sum(x) for x in history]
    push!(overall_freqs, coop_freqs)

    plt.plot(collect(1:numgens+1), coop_freqs, label="run number $i",
        c = clr[i])
end

ax = fig.gca()
ax.legend(loc=2)
ax.set_xscale("log")
ax.set_xlabel("time")
ax.set_ylabel(L"\rho_c")
ax.set_xlim([1,numgens])
fig.tight_layout()
plt.gcf()

# anim = @animate for i=1:length(history)
#     heatmap(history[i])
# end

#
# macro_lattices = []
# coop_freqs = []
#
# fig = figure()
#
# c_vals = [0.05, 0.25, 0.45, 0.65, 0.85]
# for c in c_vals
#     lattices = []
#     push!(lattices, lattice)
#     for i in 1:numgens-1
#         #println("$lattice")
#         println("$i")
#         global lattice = evolve(lattice, payoff_matrix, c, κ)
#         #imshow(lattice)
#         gcf()
#         push!(lattices, lattice)
#     end
#
#     coop_freq = [0.5]
#     for i in 1:numgens
#         push!(coop_freq, sum(lattices[i])/N^2)
#     end
#     plot(coop_freq, label="c = $c")
#     push!(macro_lattices, lattices)
#     push!(coop_freqs, coop_freq)
# end
# xscale("log")
# xlim([1,numgens])
# legend(loc=2)
# xlabel("time")
# ylabel("cooperator frequency")
# fig.savefig("test_coop_freqs.pdf")
#
# for (ci, c) in enumerate(c_vals)
#     figure()
#     anim = @animate for i=1:length(macro_lattices[ci])
#         heatmap(macro_lattices[ci][i])
#     end
#     gif(anim, "lattice_c_" * "$c" * "_fps15_take2.gif", fps = 15)
# end
