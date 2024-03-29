#!/usr/bin/env julia

## pair_ODEs.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Numerically solve the pair approximation equations
## as in Li et al. (2019).

using Revise
using DifferentialEquations
using LatticeCoop
using LaTeXStrings
using PyPlot
#pyplot()

# this commented block is just to convince me that p_cd := p_dc

# pop.neighbor_pairs = LatticeCoop.get_neighbor_pairs(pop.lattice,
#     pop.neighbors, false)
#
#
# pair_probs = fill(Float64[0], 4)
#
# ρ = [Float64[1.0 - sum(pop.lattice)/N^2], Float64[1.0*sum(pop.lattice)/N^2]]
#
# for (pi, pair) in enumerate(pop.neighbor_pairs)
#     if !pop.lattice[pair[1]] && !pop.lattice[pair[2]]
#         pair_probs[1][1] += 1.0/4/N^2
#     elseif !pop.lattice[pair[1]] && pop.lattice[pair[2]]
#         pair_probs[2][1] += 1.0/4/N^2
#     elseif pop.lattice[pair[1]] && !pop.lattice[pair[2]]
#         pair_probs[3][1] += 1.0/4/N^2
#     else
#         pair_probs[4][1] += 1.0/4/N^2
#     end
# end

function transaction_prob(
    i::Int64,
    j::Int64
    )
    # compute the likelihood of neighbors i and j paying
    # a certain number of transaction costs
    if i == 0 && j == 0
        # this can't happen: one must pay the other
        return 0.0
    elseif i == 4 && j == 4
        # nor can this: both of them can't pay
        return 0.0
    elseif i == 0 || j == 4
        # this means j must have paid i
        return binomial(3, i) * binomial(3, j-1) * 0.5^6 * 0.5
    elseif i == 4 || j == 0
        # this means i must have paid j
        return binomial(3, i-1) * binomial(3, j) * 0.5^6 * 0.5
    else
        # either x or y could have paid the other
        # so consider both possibilities
        return binomial(3, i-1) * binomial(3, j) * 0.5^6 * 0.5 +
            binomial(3, i) * binomial(3, j-1) * 0.5^6 * 0.5
    end
end

function H(
    indv_i::Bool,
    indv_j::Bool,
    i_neighbors::Tuple{Bool, Bool, Bool},
    j_neighbors::Tuple{Bool, Bool, Bool},
    game::LatticeGame
    )
    # compute the average fermi function ("energy"/update probability)
    # between two neighbors
    # by summing over their (specified) neighbors
    # and summing over all possible transaction arrangements

    # N.B.: this should return the probability that indv_j
    # will adopt indv_i's strategy
    # or, alternately, that indv_i will invade indv_j

    # this converts strategy Bool into a vector
    i_strat = [indv_i == true; indv_i == false]
    j_strat = [indv_j == true; indv_j == false]

    i_payoff = 0
    j_payoff = 0

    for (k, k_val) in enumerate(i_neighbors)
        neighbor_strat = [k_val == true; k_val == false]
        # don't worry about the cost
        # we'll deal with that later
        i_payoff += transpose(i_strat)*game.A*neighbor_strat
    end
    i_payoff += transpose(i_strat)*game.A*j_strat
    for (k, k_val) in enumerate(j_neighbors)
        neighbor_strat = [k_val == true; k_val == false]
        j_payoff += transpose(j_strat)*game.A*neighbor_strat
    end
    j_payoff += transpose(j_strat)*game.A*i_strat


    average_energy = 0

    # this is where we deal with the cost
    # basically, we iterate over each possible transaction arrangement
    # weighted by their probability
    # this is a Binomial(3, 0.5) for each individual
    # plus a coin flip to determine who paid whom
    for i_trans in 0:4
        for j_trans in 0:4
            prob = transaction_prob(i_trans, j_trans)
            payoff_exponent = -((i_payoff - game.c*i_trans) -
                (j_payoff - game.c*j_trans))/game.κ
            average_energy += prob/(1.0+exp(payoff_exponent))
        end
    end
    return average_energy
end

function p(
    strat1::Bool,
    strat2::Bool,
    pair_probs::Array{Float64, 1}
    )
    # this returns the p_xy pair probabilities within the sums
    return pair_probs[2*strat1+strat2+1]
end

function dot_pair_probs(
    d_pair_probs::Array{Float64, 1},
    pair_probs::Array{Float64, 1},
    game::LatticeGame,
    t::Float64
    )
    # println("$t, $pair_probs")
    # this function is a mouthful
    # it implements equations 6 and 7 from Li et al. (2019)
    # pair_probs indices are as follows: dd, dc, cd, cc

    # compute the frequencies of c and d
    (ρ_c, ρ_d) = (pair_probs[4] + pair_probs[2], pair_probs[3] + pair_probs[1])

    # multiplicative term preceding each sum
    # the 2 comes from symmetry
    # the ρ terms are normalizing constants
    # note that this prefactor does not affect equilibrium dynamics
    prefactor = 2*pair_probs[2]/(ρ_c^3*ρ_d^3)

    # get a list of possible neighbor configurations
    # there will be 3 of each and they can be true or false
    i_neighbor_configs = vec(collect(Base.Iterators.product(fill([true, false], 3)...)))
    j_neighbor_configs = vec(collect(Base.Iterators.product(fill([true, false], 3)...)))

    # compute the change in p_cc
    outer_sum = 0
    # first term
    for (ii, i_neighbors) in enumerate(i_neighbor_configs)
        # count the number of neighbor cooperators
        n_c = sum(i_neighbors)
        # start computing the inner sum
        inner_sum = 0
        # sum over all possible values of u, v, and w
        # we need p_cu, p_cv, and p_cw for this
        for (ji, j_neighbors) in enumerate(j_neighbor_configs)
            inner_sum += prod([p(true, j_neighbors[k], pair_probs) for k in 1:3]) *
                H(false, true, i_neighbors, j_neighbors, game)
                # the above should give H[P_c(u,v,w) → P_d(x,y,z)]
                # or, in other words, the probability that
                # j, a cooperator with neighbors (u, v, w)
                # will adopt the strategy of
                # i, a defector with neighbors (x, y, z)
        end
        # compute the outer sum and multiply terms by the inner sum
        # we need p_dx, p_dy, and p_dz
        outer_sum += (n_c + 1) *
            prod([p(false, i_neighbors[k], pair_probs) for k in 1:3]) *
            inner_sum
    end

    # second term
    for (ii, i_neighbors) in enumerate(i_neighbor_configs)
        # count the number of neighbor cooperators
        n_c = sum(i_neighbors)
        # start computing the inner sum
        inner_sum = 0
        # sum over all possible values of u, v, and w
        # we need p_du, p_dv, and p_dw for this
        for (ji, j_neighbors) in enumerate(j_neighbor_configs)
            inner_sum += prod([p(false, j_neighbors[k], pair_probs) for k in 1:3]) *
                H(true, false, i_neighbors, j_neighbors, game)
                # the above should give H[P_d(u,v,w) → P_c(x,y,z)]
                # or, in other words, the probability that
                # j, a defector with neighbors (u, v, w)
                # will adopt the strategy of
                # i, a cooperator with neighbors (x, y, z)
        end
        # compute the outer sum and multiply terms by the inner sum
        # we need p_cx, p_cy, and p_cz
        outer_sum -= n_c *
            prod([p(true, i_neighbors[k], pair_probs) for k in 1:3]) *
            inner_sum
    end

    # this gives d/dt(p_cc)
    d_pair_probs[4] = prefactor*outer_sum


    # compute the change in p_cd
    outer_sum = 0
    # first term
    for (ii, i_neighbors) in enumerate(i_neighbor_configs)
        # count the number of neighbor cooperators
        n_c = sum(i_neighbors)
        # start computing the inner sum
        inner_sum = 0
        # sum over all possible values of u, v, and w
        # we need p_cu, p_cv, and p_cw for this
        for (ji, j_neighbors) in enumerate(j_neighbor_configs)
            inner_sum += prod([p(true, j_neighbors[k], pair_probs) for k in 1:3]) *
                H(false, true, i_neighbors, j_neighbors, game)
                # the above should give H[P_c(u,v,w) → P_d(x,y,z)]
                # or, in other words, the probability that
                # j, a cooperator with neighbors (u, v, w)
                # will adopt the strategy of
                # i, a defector with neighbors (x, y, z)
        end
        # compute the outer sum and multiply by the inner sum
        # we need p_dx, p_dy, and p_dz for this
        outer_sum += (1 - n_c) *
            prod([p(false, i_neighbors[k], pair_probs) for k in 1:3]) *
            inner_sum
    end
    # second term
    for (ii, i_neighbors) in enumerate(i_neighbor_configs)
        # count the number of neighbor cooperators
        n_c = sum(i_neighbors)
        # start computing the inner sum
        inner_sum = 0
        # sum over all possible values of u, v, and w
        # we need p_du, p_dv, and p_dw for this
        for (ji, j_neighbors) in enumerate(j_neighbor_configs)
            inner_sum += prod([p(false, j_neighbors[k], pair_probs) for k in 1:3]) *
                H(true, false, i_neighbors, j_neighbors, game)
                # the above should give H[P_d(u,v,w) → P_c(x,y,z)]
                # or, in other words, the probability that
                # j, a defector with neighbors (u, v, w)
                # will adopt the strategy of
                # i, a cooperator with neighbors (x, y, z)
        end
        # compute the outer sum and multiply terms by the inner sum
        # we need p_cx, p_cy, and p_cz
        outer_sum -= (2 - n_c) *
            prod([p(true, i_neighbors[k], pair_probs) for k in 1:3]) *
            inner_sum
    end
    # that should do it

    # by symmetry, the following two terms are the same
    d_pair_probs[3] = prefactor*outer_sum # d/dt(p_cd)
    d_pair_probs[2] = prefactor*outer_sum # d/dt(p_dc)

    # the terms should all sum to zero, so this gives d/dt(p_dd)
    d_pair_probs[1] = -sum(d_pair_probs[2:4])
end

function dot_pair_probs_revised(
    d_pair_probs::Array{Float64, 1},
    pair_probs::Array{Float64, 1},
    game::LatticeGame,
    t::Float64
    )
    # this function is a mouthful
    # it implements equations 6 and 7 from Li et al. (2019)
    # pair_probs indices are as follows: dd, dc, cd, cc
    # here i considered the possibility that some of the indices
    # have accidentally been switched
    # (this doesn't appear to be the case)

    # compute the frequencies of c and d
    (ρ_c, ρ_d) = (pair_probs[4] + pair_probs[2], pair_probs[3] + pair_probs[1])

    # multiplicative term preceding each sum
    # the 2 comes from symmetry
    # the ρ terms are normalizing constants
    # note that this prefactor does not affect equilibrium dynamics
    prefactor = 2*pair_probs[2]/(ρ_c^3*ρ_d^3)

    # get a list of possible neighbor configurations
    # there will be 3 of each and they can be true or false
    i_neighbor_configs = vec(collect(Base.Iterators.product(fill([true, false], 3)...)))
    j_neighbor_configs = vec(collect(Base.Iterators.product(fill([true, false], 3)...)))

    # compute the change in p_cc
    outer_sum = 0
    # first term
    for (ii, i_neighbors) in enumerate(i_neighbor_configs)
        # count the number of neighbor cooperators
        # start computing the inner sum
        inner_sum = 0
        # sum over all possible values of u, v, and w
        # we need p_cu, p_cv, and p_cw for this
        for (ji, j_neighbors) in enumerate(j_neighbor_configs)
            n_c = sum(j_neighbors)
            inner_sum += (n_c + 1) *
                prod([p(true, j_neighbors[k], pair_probs) for k in 1:3]) *
                H(false, true, i_neighbors, j_neighbors, game)
                # the above should give H[P_c(u,v,w) → P_d(x,y,z)]
                # or, in other words, the probability that
                # j, a cooperator with neighbors (u, v, w)
                # will adopt the strategy of
                # i, a defector with neighbors (x, y, z)
        end
        # compute the outer sum and multiply terms by the inner sum
        # we need p_dx, p_dy, and p_dz
        outer_sum +=
            prod([p(false, i_neighbors[k], pair_probs) for k in 1:3]) *
            inner_sum
    end

    # second term
    for (ii, i_neighbors) in enumerate(i_neighbor_configs)
        # count the number of neighbor cooperators
        # start computing the inner sum
        inner_sum = 0
        # sum over all possible values of u, v, and w
        # we need p_du, p_dv, and p_dw for this
        for (ji, j_neighbors) in enumerate(j_neighbor_configs)
            n_c = sum(j_neighbors)
            inner_sum += n_c *
                prod([p(false, j_neighbors[k], pair_probs) for k in 1:3]) *
                H(true, false, j_neighbors, i_neighbors, game)
                # the above should give H[P_d(u,v,w) → P_c(x,y,z)]
                # or, in other words, the probability that
                # j, a defector with neighbors (u, v, w)
                # will adopt the strategy of
                # i, a cooperator with neighbors (x, y, z)
        end
        # compute the outer sum and multiply terms by the inner sum
        # we need p_cx, p_cy, and p_cz
        outer_sum -=
            prod([p(true, i_neighbors[k], pair_probs) for k in 1:3]) *
            inner_sum
    end

    # this gives d/dt(p_cc)
    d_pair_probs[4] = prefactor*outer_sum


    # compute the change in p_cd
    outer_sum = 0
    # first term
    for (ii, i_neighbors) in enumerate(i_neighbor_configs)
        # count the number of neighbor cooperators
        # start computing the inner sum
        inner_sum = 0
        # sum over all possible values of u, v, and w
        # we need p_cu, p_cv, and p_cw for this
        for (ji, j_neighbors) in enumerate(j_neighbor_configs)
            n_c = sum(j_neighbors)
            inner_sum += (1 - n_c) *
                prod([p(true, j_neighbors[k], pair_probs) for k in 1:3]) *
                H(false, true, i_neighbors, j_neighbors, game)
                # the above should give H[P_c(u,v,w) → P_d(x,y,z)]
                # or, in other words, the probability that
                # j, a cooperator with neighbors (u, v, w)
                # will adopt the strategy of
                # i, a defector with neighbors (x, y, z)
        end
        # compute the outer sum and multiply by the inner sum
        # we need p_dx, p_dy, and p_dz for this
        outer_sum +=
            prod([p(false, i_neighbors[k], pair_probs) for k in 1:3]) *
            inner_sum
    end
    # second term
    for (ii, i_neighbors) in enumerate(i_neighbor_configs)
        # count the number of neighbor cooperators
        # start computing the inner sum
        inner_sum = 0
        # sum over all possible values of u, v, and w
        # we need p_du, p_dv, and p_dw for this
        for (ji, j_neighbors) in enumerate(j_neighbor_configs)
            n_c = sum(j_neighbors)
            inner_sum += (2 - n_c) *
                prod([p(false, j_neighbors[k], pair_probs) for k in 1:3]) *
                H(true, false, i_neighbors, j_neighbors, game)
                # the above should give H[P_d(u,v,w) → P_c(x,y,z)]
                # or, in other words, the probability that
                # j, a defector with neighbors (u, v, w)
                # will adopt the strategy of
                # i, a cooperator with neighbors (x, y, z)
        end
        # compute the outer sum and multiply terms by the inner sum
        # we need p_cx, p_cy, and p_cz
        outer_sum -=
            prod([p(true, i_neighbors[k], pair_probs) for k in 1:3]) *
            inner_sum
    end
    # that should do it

    # by symmetry, the following two terms are the same
    d_pair_probs[3] = prefactor*outer_sum # d/dt(p_cd)
    d_pair_probs[2] = prefactor*outer_sum # d/dt(p_dc)

    # the terms should all sum to zero, so this gives d/dt(p_dd)
    d_pair_probs[1] = -sum(d_pair_probs[2:4])
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
