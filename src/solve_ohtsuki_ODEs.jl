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
    d_x::Array{Float64, 1}, # \dot{x}_i
    x::Array{Float64, 1}, # current x_i
    params::Array{Any, 1},
    t::Float64, # time
    # k::Int64, # degree of network
    # w::Float64, # strength of selection
    # a::Array{Float64, 2}, # game
    # update_rule::String="birth_death" # update rule
    )

    k, w, a, update_rule = params

    num_strats = size(a)[1] # get the number of strategies
    f = zeros(num_strats)
    g = zeros(num_strats)
    b = zeros(num_strats, num_strats)


    # populate the f, g, and b arrays

    # f_i = \sum_j x_j a_{ij}
    [f[i] = sum([x[j] * a[i,j] for j in 1:num_strats]) for i in 1:num_strats]

    if update_rule == "birth_death" || update_rule == "bd"
        # b_{ij} = (a_{ii} + a_{ij} - a_{ji} - a_{jj})/(k-2)
        [[b[i,j] = a[i,i] + a[i,j] -
            a[j,i] - a[j,j] for j in 1:num_strats] for i in 1:num_strats]
        b /= (k-2)
    end

    # g_i = \sum_k x_j b_{ij}
    [g[i] = sum([x[j] * b[i,j] for j in 1:num_strats]) for i in 1:num_strats]

    # phi = \sum_i x_i f_i = \sum_{i,j} x_i x_j a_{ij}
    phi = sum([sum([x[i]*x[j]*a[i,j] for j in 1:num_strats]) for i in 1:num_strats])

    # \dot{x}_i = x_i (f_i + g_i - phi)
    # note that we need to edit the existing d_x array here rather than
    # create a new one!
    [d_x[i] = w*(k-2)^2/(k-1) * x[i] * (f[i] + g[i] - phi) for i in 1:num_strats]
    println("$k, $w, $a, $f, $g, $b, $phi, $x, $d_x")

end

k = 5
w = 0.1
a = [5 0; 8 1]
max_time = 100.0

num_strats = size(a)[1]

# initial pair probabilities
freqs_0 = fill(1.0/num_strats, num_strats)
#freqs_0 = [0.1, 0.9]
tspan = (0.0, max_time)
# solve the ODE
params = [k, w, a, "bd"]
prob = ODEProblem(dot_x_i, freqs_0, tspan, params)
sol = solve(prob)


# function solve_ODEs(
#     k::Float64,
#     w::Float64,
#     a::Array{Float64, 2},
#     update_rule::String,
#     max_time::Float64,
#     make_plot::Bool = false
#     )
#
#     num_strats = size(a)[1]
#
#     # initial pair probabilities
#     freqs_0 = fill(1.0/num_strats, num_strats)
#     tspan = (0.0, max_time)
#     # solve the ODE
#     prob = ODEProblem(dot_x_i, freqs_0, tspan, k, w, a)
#     sol = solve(prob)
#     if make_plot
#         fig = plt.figure()
#         labels = [L"p_{dd}", L"p_{dc}", L"p_{cd}", L"p_{cc}"]
#         for y in reverse(1:4)
#             plt.plot(sol.t, [sol.u[x][y] for x in 1:(length(sol.t))], label=labels[y])
#         end
#         legend(loc=4)
#         ax=plt.gca()
#         ax.set_xlabel("time")
#         ax.set_ylabel("frequency")
#         ax.set_ylim([0.0,1.0])
#         plt.title("b = $b, c = $c, κ = $κ")
#         plt.tight_layout()
#         display(fig)
#     end
#     # this has several attributes, including the time sol.t
#     # and the actual solution sol.u
#     return sol
# end
#
# max_time = 20.0
# make_plots = true
#
# b_step = 0.01
# c_step = 0.1
#
# b_vals = collect(1.0:b_step:1.06)
# c_vals = collect(0.0:c_step:1.0)
# κ = 0.1
# ρ_c_vals = zeros(length(b_vals), length(c_vals))
# for (bi, b) in enumerate(b_vals)
#     for (ci, c) in enumerate(c_vals)
#         println("b = $b, c = $c")
#         sol = solve_ODEs(b, c, κ, max_time, make_plots)
#         ρ_c_vals[bi, ci] = sum(sol.u[end][3:4])
#     end
# end
#
# fig = plt.figure()
# color_list = ["red", "blue", "yellow", "green", "purple", "gray", "pink"]
# for (bi, b) in enumerate(b_vals)
#     plt.plot(c_vals, ρ_c_vals[bi,:], label="b = $b", c = color_list[bi])
# end
# ax=plt.gca()
# ax.set_xlabel("c")
# ax.set_ylabel(L"\rho_c")
# ax.set_ylim([0.0,0.7])
# ax.set_xlim([0.0,1.0])
# plt.legend(loc=3)
# plt.tight_layout()
# display(fig)
