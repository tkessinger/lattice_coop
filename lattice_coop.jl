#!/usr/bin/env julia

## lattice_coop.jl
##
## Author: Taylor Kessinger <taylor.kessinger@uky.edu>
## First stab at duplicating Li et al. (2019).

using Random
using PyPlot
using Revise
using LaTeXStrings

# IDEA: I can't imagine that computing the energy gap between every interacting
# pair is necessary or efficient.
# This step could probably be sped up by computing the possible energy values
# in advance, then calling them on the fly.
# This would be motivated by the fact that there are really only a few
# possible values the payoff can take.
# We would need to track the strategy of each individual's neighbors
# and in how many games they and their neighbors paid the transaction cost.
# Then, simply look up the energy values in a table.
#
# Update: It turns out "a few" is on the order of 400.
# Maybe this idea is doomed to failure.
# Must profile.

# Organizational to-dos:
#   Make this a module. Probably have a "Lattice" type and a "Game" type.
#   (This would cut down on parameter passing.)
#   Fix that annoying global lattice update. It's bad form.
#   Lattice type might have neighbors built in.
#   This might be unnecessary, though, and would make iterating harder.
#   Figure out a smarter way to save histories,
#   e.g., a "History" type or attribute with BitArrays.
#   That should be doable since BitArray is so space efficient.
#   A 100x100 BitArray should in principle use only

struct LatticeGame
    # new type for "game" parameters
    # this should minimize the amount of passing around that needs to be done

    b::Float64 # benefit to defecting
    c::Float64 # transaction cost
    κ::Float64 # temperature
    A::Array{Float64,2} # game structure

    # constructors
    LatticeGame(b,c,κ) = new(b,c,κ,Float64[1 0; b 0])
    # [1 0; b 0] is the form for the "simplified" prisoner's dilemma game
end

# mutable struct Lattice
#     # new type for lattices
#     # this should make it easier to pass around neighbors, etc.
#
#     N::Int64 # size of lattice
#     lattice::BitArray # false are defectors, true are cooperators
#     neighbors::Dict{CartesianIndex,Array{CartesianIndex,1}} # list of neighbors for each individual
#     game::LatticeGame # game parameters
#
# end

function assign_transactions(lattice::BitArray)
    # decide on an ordering for individuals
    # individuals with higher order will always initiate transactions
    # against individuals with lower order
    # (and hence pay the cost)
    (width, height) = size(lattice)
    ordering = reshape(randperm(height*width), height, width)
    if verbose println("ordering = $ordering") end
    return ordering
end

function get_strategy(bool::Bool)
    # converts from true-false (coop-defect)
    # into vector form
    return [bool==true; bool==false]
end

function strategy_string(strat::Array{Bool})
    if strat == [true; false]
        return "cooperator"
    elseif strat == [false; true]
        return "defector"
    else
        println("something has gone wrong")
    end
end

function get_payoffs(lattice::BitArray, neighbors::Dict{CartesianIndex, Array{CartesianIndex,1}},
    ordering::Array{Int64,2}, game::LatticeGame)
    # determine the payoff for every individual in the lattice
    payoffs = zeros(size(lattice))
    for (i, indv) in enumerate(CartesianIndices(lattice))
        # iterate over each individual
        tmp_payoffs = 0
        lattice_val = lattice[indv]
        for (ni, neighbor) in enumerate(neighbors[indv])
            # iterate over all their neighbors
            # get everyone's strategy
            neighbor_val = lattice[neighbor]
            indv_strat = get_strategy(lattice_val)
            neighbor_strat = get_strategy(neighbor_val)

            # decide who has to pay the transaction cost
            pay_cost = ordering[indv] > ordering[neighbor]
            tmp_payoff = transpose(indv_strat)*game.A*neighbor_strat - game.c*(pay_cost)

            if verbose
                println("individual $(indv.I) has strategy $(strategy_string(indv_strat)) and ordering $(ordering[indv])")
                println("neighbor $(neighbor.I) has strategy $(strategy_string(neighbor_strat)) and ordering $(ordering[neighbor])")
                println("individual should pay transaction cost: $pay_cost")
                println("payoff is $tmp_payoff")
            end

            tmp_payoffs += tmp_payoff
        end
        payoffs[indv] = tmp_payoffs
        if verbose println("individual $(indv.I) had net payoff $tmp_payoffs") end
    end
    return payoffs
end

function get_neighbors(lattice::BitArray, wraparound::Bool=true)
    # obtain the neighbors for every individual in the lattice
    # this only needs to be run once, at the beginning
    # returns a Dict with a CartesianIndex() key corresponding to the individual
    # value is a list of CartesianIndices() corresponding to neighbors

    # unnecessary but allows for non-square lattices
    (width, height) = size(lattice)

    neighbors = Dict{CartesianIndex,Array{CartesianIndex,1}}()
    for (i, indv) in enumerate(CartesianIndices(lattice))
        # this loop just controls for the possibility that an individual
        # is along an edge
        # if it is, and wraparound is set, it handles that accordingly
        # otherwise, it grants fewer neighbors
        tmp_neighbors = CartesianIndex[]
        (x, y) = (indv[1], indv[2])
        if x != 1
            push!(tmp_neighbors, CartesianIndex(x-1, y))
        end
        if y != 1
            push!(tmp_neighbors, CartesianIndex(x, y-1))
        end
        if x != width
            push!(tmp_neighbors, CartesianIndex(x+1, y))
        end
        if y != height
            push!(tmp_neighbors, CartesianIndex(x, y+1))
        end

        # i am absolutely certain there is a less ugly way to do this
        if width > 2
            if x == 1 && wraparound
                push!(tmp_neighbors, CartesianIndex(width, y))
            end
            if x == width && wraparound
                push!(tmp_neighbors, CartesianIndex(1, y))
            end
        end
        if height > 2
            if y == 1 && wraparound
                push!(tmp_neighbors, CartesianIndex(x, height))
            end
            if y == height && wraparound
                push!(tmp_neighbors, CartesianIndex(x, 1))
            end
        end

        neighbors[indv] = tmp_neighbors
    end
    return neighbors
end

function energy_function(indv::CartesianIndex, neighbor::CartesianIndex,
    lattice::BitArray, payoffs::Array{Float64}, game::LatticeGame)
    # return the "energy" difference between individuals
    # this sets the probability that one of them flips to the other's strategy
    payoff_diff = payoffs[indv] - payoffs[neighbor]
    energy = 1.0/(1.0+exp(payoff_diff/game.κ))
    if !(energy in observed_energy_values)
        push!(observed_energy_values, energy)
    end
    return energy
end

function evolve(lattice::BitArray, game::LatticeGame)
    # evolves the lattice one generation
    # determines an ordering, then determines payoffs
    # (including transaction costs)
    # then updates the lattice as individuals switch strategies

    ordering = assign_transactions(lattice)
    payoffs = get_payoffs(lattice, neighbors, ordering, game)
    new_lattice = BitArray(zeros(size(lattice)))

    # unnecessary but allows for a non-square lattice
    height, width = size(lattice)

    update = rand(height, width)

    if verbose println("update array looks like $update") end

    for (i, indv) in enumerate(CartesianIndices(lattice))
        # pick a random neighbor for each individual and compare energies
        # the energy sets the probability of a state flip
        neighbor = rand(neighbors[indv])
        if verbose
            println("individual $(indv.I) chose random neighbor $(neighbor.I)")
            println("individual $(indv.I) had payoff $(payoffs[indv])")
            println("neighbor $(neighbor.I) had payoff $(payoffs[neighbor])")
            println("the energy function is $(energy_function(indv, neighbor, lattice, payoffs, game))")
            println("the update function is $(update[indv])")
        end
        if update[indv] < energy_function(indv, neighbor, lattice, payoffs, game)
            new_lattice[indv] = lattice[neighbor]
            if verbose println("so individual $(indv.I) updates their strategy") end
        else
            new_lattice[indv] = lattice[indv]
            if verbose println("so individual $(indv.I) keeps their strategy") end
        end
    end
    lattice = new_lattice
    return lattice
end

b = 1.02 # benefit to defecting
c = 0.26 # cost for initiating transaction
κ = 0.1 # temperature parameter

verbose = false # turn this on for detailed error tracking

histories = []
overall_freqs = []

fig = figure()

c_list = [0.05, 0.25, 0.45, 0.65, 0.85]
for (ci, c) in enumerate(c_list)

    println("evolving with transaction cost $c")

    numgens = 1000
    N = 100 # size (dimension) of lattice
    game = LatticeGame(b, c, κ)


    # initialize a lattice with zeros
    # individuals will be specified by a CartestianIndex()
    # of coordinates on the lattice
    lattice = BitArray(zeros(N,N))
    # randomly flip half of them to ones
    rand_init = rand(CartesianIndices(lattice), floor(Int64, N^2/2))
    [lattice[x] = true for x in rand_init]

    # determine every individual's neighbors
    # i am sure there is a smarter way to do this
    neighbors = get_neighbors(lattice)

    history = []
    push!(history, lattice)
    for i in 1:numgens
        lattice = evolve(lattice, game)
        push!(history, lattice)
    end

    push!(histories, history)

    coop_freqs = 1.0/N^2*[sum(x) for x in history]
    push!(overall_freqs, coop_freqs)

    plot(coop_freqs, label="c = $c")
end
ax = fig.gca()
ax.legend(loc=2)
ax.set_xscale("log")
ax.set_xlabel("time")
ax.set_ylabel(L"\rho_c")
fig.tight_layout()
gcf()

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
