#!/usr/bin/env julia

## solve_ohtsuki_ODEs.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Numerically solve the pair approximation equations
## as in Ohtsuki et al. (2006).

using Revise
using DifferentialEquations
using LatticeCoop
using LaTeXStrings
using PyPlot
#pyplot()

function dot_x_i(
    d_x::Array{Float64, 1}, #\dot{x}_i
    x::Array{Float64, 1}, # current x_i
    t::Float64, # time
    k::Int64, # degree of network
    w::Float64, # strength of selection
    a::Array{Float64, 2}, # game
    update_rule::String="birth_death" # update rule
    )
    num_strats = size(game)[1]
    f = zeros(num_strats)
    g = zeros(num_strats)
    b = zeros(num_strats, num_strats)

    [f[i] = sum([x[j] * a[i,j] for j in 1:num_strats]) for i in 1:num_strats]
    if update_rule == "birth_death" || update_rule == "bd"
        [[b[i,j] = a[i,i] + a[i,j] -
            a[j,i] - a[j,j] for j in 1:num_strats] for i in 1:num_strats]
        b /= (k-2)
    end
    [g[i] = sum([x[j] * b[i,j] for j in 1:num_strats]) for i in 1:num_strats]
    phi = sum([sum([x[i]*x[j]*a[i,j] for j in 1:num_strats]) for i in 1:num_strats])

end

function solve_ODEs(
    b::Float64,
    c::Float64,
    κ::Float64,
    max_time::Float64,
    make_plot::Bool = false
    )

    # initialize a LatticeGame
    game = LatticeGame(b, c, κ)
    # initial pair probabilities
    pair_probs_0 = [0.25, 0.25, 0.25, 0.25]
    tspan = (0.0, max_time)
    # solve the ODE
    prob = ODEProblem(dot_pair_probs_revised, pair_probs_0, tspan, game)
    sol = solve(prob)
    if make_plot
        fig = plt.figure()
        labels = [L"p_{dd}", L"p_{dc}", L"p_{cd}", L"p_{cc}"]
        for y in reverse(1:4)
            plt.plot(sol.t, [sol.u[x][y] for x in 1:(length(sol.t))], label=labels[y])
        end
        legend(loc=4)
        ax=plt.gca()
        ax.set_xlabel("time")
        ax.set_ylabel("frequency")
        ax.set_ylim([0.0,1.0])
        plt.title("b = $b, c = $c, κ = $κ")
        plt.tight_layout()
        display(fig)
    end
    # this has several attributes, including the time sol.t
    # and the actual solution sol.u
    return sol
end

max_time = 20.0
make_plots = true

b_step = 0.01
c_step = 0.1

b_vals = collect(1.0:b_step:1.06)
c_vals = collect(0.0:c_step:1.0)
κ = 0.1
ρ_c_vals = zeros(length(b_vals), length(c_vals))
for (bi, b) in enumerate(b_vals)
    for (ci, c) in enumerate(c_vals)
        println("b = $b, c = $c")
        sol = solve_ODEs(b, c, κ, max_time, make_plots)
        ρ_c_vals[bi, ci] = sum(sol.u[end][3:4])
    end
end

fig = plt.figure()
color_list = ["red", "blue", "yellow", "green", "purple", "gray", "pink"]
for (bi, b) in enumerate(b_vals)
    plt.plot(c_vals, ρ_c_vals[bi,:], label="b = $b", c = color_list[bi])
end
ax=plt.gca()
ax.set_xlabel("c")
ax.set_ylabel(L"\rho_c")
ax.set_ylim([0.0,0.7])
ax.set_xlim([0.0,1.0])
plt.legend(loc=3)
plt.tight_layout()
display(fig)
