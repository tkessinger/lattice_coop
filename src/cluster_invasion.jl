#!/usr/bin/env julia

## cluster_invasion.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot fitness and update functions and their dependence on c.
## This should give us some intuition as to why invasion
## becomes easier or harder at different c values.

using Revise
using LatticeCoop
using PyPlot

c = 0.26
κ = 0.1

function H(
    i_payoff::Float64,
    j_payoff::Float64,
    i_trans::Float64,
    j_trans::Float64,
    game::LatticeGame,
    )
    return ((i_payoff - game.c*i_trans) -
        (j_payoff - game.c*j_trans))/game.κ
end

function stats_H(
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
    # returns both the average of H
    # and the variance of H
    # as well as the average and variance of i and j's fitness

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

    average_fitness_i = 0
    var_fitness_i = 0
    average_fitness_j = 0
    var_fitness_j = 0

    average_energy = 0
    var_energy = 0

    # this is where we deal with the cost
    # basically, we iterate over each possible transaction arrangement
    # weighted by their probability
    # this is a Binomial(3, 0.5) for each individual
    # plus a coin flip to determine who paid whom
    for i_trans in 0:4
        for j_trans in 0:4
            prob = transaction_prob(i_trans, j_trans)
            fitness_i = i_payoff - game.c*i_trans
            fitness_j = j_payoff - game.c*j_trans
            payoff_exponent = (fitness_i -
                fitness_j)/game.κ
            average_fitness_i += prob*fitness_i
            average_fitness_j += prob*fitness_j
            var_fitness_i += prob*fitness_i^2
            var_fitness_j += prob*fitness_j^2

            average_energy += prob/(1.0+exp(payoff_exponent))
            var_energy += prob/(1.0+exp(payoff_exponent))^2
        end
    end
    var_fitness_i -= average_fitness_i^2
    var_fitness_j -= average_fitness_j^2
    var_energy -= average_energy^2
    return [average_fitness_i, var_fitness_i,
    average_fitness_j, var_fitness_j,
    average_energy, var_energy]
end

function strategy_string(
    strat::Bool,
    full_name::Bool=false
    )
    if strat==false && full_name return "defector"
    elseif strat==false return "d"
    elseif strat==true && full_name return "cooperator"
    elseif strat==true return "c" end
end

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

n_c_i = 1
n_c_j = 1

c_step = 0.01
b_step = 0.02

c_list = 0.0:c_step:1.0
b_list = 1.0:b_step:1.06


i_strat = true
j_strat = true

H_array = zeros(4,4,length(c_list), length(b_list), 6)
# last index:
# 1. E(F_i)
# 2. Var(F_i)
# 3. E(F_j)
# 4. Var(F_j)
# 5. E(H)
# 6. Var(H)

κ = 0.1

n_c_i_list = 0:3
n_c_j_list = 0:3
for n_c_i in n_c_i_list
    for n_c_j in n_c_j_list
        i_neighbors = (0 < n_c_i, 1 < n_c_i, 2 < n_c_i)
        j_neighbors = (0 < n_c_j, 1 < n_c_j, 2 < n_c_j)
        for (ci, c_val) in enumerate(c_list)
            for (bi, b_val) in enumerate(b_list)
                game = LatticeGame(b_val, c_val, κ)
                H_array[n_c_i+1,n_c_j+1,ci,bi,:] =
                    stats_H(i_strat, j_strat, i_neighbors, j_neighbors, game)
            end
        end
    end
end

labels = [L"\bar{F_i}",
    L"\mathrm{Var}(F_i)",
    L"\bar{F_j}",
    L"\mathrm{Var}(F_i)",
    L"\bar{H}",
    L"\mathrm{Var}(H)"]
save_strings = ["avg_Fi",
    "var_Fi",
    "avg_Fj",
    "var_Fj",
    "avg_H",
    "var_H"]

for k in 1:6
    fig, axs = plt.subplots(4,4,figsize=(10,10), sharex="all", sharey="all")
    for n_c_i in n_c_i_list
        for n_c_j in n_c_j_list
            tmp_ax = axs[n_c_i+1,n_c_j+1]
            for (bi, b_val) in enumerate(b_list)
                #tmp_ax.vlines(0.26,0,.15, ls="dotted")
                tmp_ax.plot(c_list, H_array[n_c_i+1,n_c_j+1,:,bi,k],label="b = $b_val")
            end
            tmp_ax.set_title(L"n_{c_i} = "*string(n_c_i) * L", n_{c_j} = "*string(n_c_j))
            #tmp_ax.set_xlabel(L"c")
            #tmp_ax.set_ylabel(labels[k])
            #tmp_ax.legend(loc=2)
        end
    end
    #for ax in axs.flat()
    #    ax.set(xlabel = L"c", ylabel = labels[k])
    #end
    [ax.set(xlabel=L"c") for ax in axs[4,:]]
    [ax.set(ylabel=labels[k]) for ax in axs[:,1]]
    axs[4,1].legend(loc=2)
    fig.suptitle("i = $(strategy_string(i_strat, true)), j = $(strategy_string(j_strat, true))")
    plt.tight_layout()
    fig.subplots_adjust(top=0.93)
    display(fig)
    plt.savefig("figures/H_stats_i_" * strategy_string(i_strat) * "_j_" *
        strategy_string(j_strat) * "_" * save_strings[k] * ".pdf")
end
# for ax in axs.flat:
#     ax.set(xlabel="x-label", ylabel='y-label')
# end
# for ax in axs.flat:
#     ax.label_outer()
# end
