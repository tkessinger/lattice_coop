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
    f = zeros(Float64, num_strats)
    g = zeros(Float64, num_strats)
    b = zeros(Float64, num_strats, num_strats)


    # populate the f, g, and b arrays

    # f_i = \sum_j x_j a_{ij}
    [f[i] = sum([x[j] * a[i,j] for j in 1:num_strats]) for i in 1:num_strats]

    if update_rule == "birth_death" || update_rule == "bd"
        # b_{ij} = (a_{ii} + a_{ij} - a_{ji} - a_{jj})/(k-2)
        [[b[i,j] = a[i,i] + a[i,j] -
            a[j,i] - a[j,j] for j in 1:num_strats] for i in 1:num_strats]
        b /= (k-2)
        prefactor = w*(k-2)^2/(k-1)
    elseif update_rule == "death_birth" || update_rule == "db"
        [[b[i,j] = (k+1)*a[i,i] + a[i,j] -
            a[j,i] - (k+1)*a[j,j] for j in 1:num_strats] for i in 1:num_strats]
        b /= (k+1)*(k-2)
        prefactor = w*(k+1)*(k-2)^2/(k*(k-1))
    elseif update_rule == "imitation" || update_rule == "im"
        [[b[i,j] = (k+3)*a[i,i] + 3*a[i,j] -
            3*a[j,i] - (k+3)*a[j,j] for j in 1:num_strats] for i in 1:num_strats]
        b /= (k+3)*(k-2)
        prefactor = w*k*(k+3)*(k-2)^2/((k+1)^2*(k-1))
    end

    # g_i = \sum_k x_j b_{ij}
    [g[i] = sum([x[j] * b[i,j] for j in 1:num_strats]) for i in 1:num_strats]

    # phi = \sum_i x_i f_i = \sum_{i,j} x_i x_j a_{ij}
    phi = sum([sum([x[i]*x[j]*a[i,j] for j in 1:num_strats]) for i in 1:num_strats])

    # \dot{x}_i = x_i (f_i + g_i - phi)
    # note that we need to edit the existing d_x array here rather than
    # create a new one!
    [d_x[i] = prefactor * x[i] * (f[i] + g[i] - phi) for i in 1:num_strats]
    println("$k, $w, $a, $f, $g, $b, $phi, $x, $d_x")

end


function solve_ODEs(
    k::Union{Int64, Float64},
    w::Union{Int64, Float64},
    a::Array{Union{Int64, Float64}, 2},
    update_rule::String,
    max_time::Union{Int64, Float64},
    make_plot::Bool = false
    )

    num_strats = size(a)[1]

    # initial frequencies
    freqs_0 = fill(1.0/num_strats, num_strats)
    tspan = (0.0, max_time)
    # solve the ODE
    params = [k, w, a, "bd"]
    prob = ODEProblem(dot_x_i, freqs_0, tspan, params)
    sol = solve(prob)
    if make_plot
        fig = plt.figure()
        labels = [L"x_c", L"x_d"]
        for y in 1:2
            plt.plot(sol.t, [sol.u[x][y] for x in 1:(length(sol.t))], label=labels[y])
        end
        legend(loc=4)
        ax=plt.gca()
        ax.set_xlabel("time")
        ax.set_ylabel("frequency")
        ax.set_ylim([0.0,1.0])
        plt.title("k = $k")#, c = $c, κ = $κ")
        plt.tight_layout()
        display(fig)
    end
    # this has several attributes, including the time sol.t
    # and the actual solution sol.u
    return sol
end


k = 3
w = 0.1
R, S, T, P = 5, 0, 8, 1
a = Float64[R S ; T P]
update_rule = "db"
max_time = 10000.0

#w_range = Float64.(10 .^(range(-3,stop=0,length=10)))
#for (wi, w) in enumerate(w_range)
    #print(w)
    solve_ODEs(k, w, a, update_rule, max_time, true)
#end
