#!/usr/bin/env julia

## check_array.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Check to make sure the different methods of determining transaction initiation
## are actually equivalent.
## (Spoiler: They aren't.
## The "ordering" method creates a uniform distribution of transactions.
## The "actually consider each pair and make one the initiator"
## doesn't--it's a binomial distribution with p=0.5.)

using Revise
using PyPlot

import LatticeCoop

N = 20
numruns = 100

# initialize an "empty" population
pop = LatticeCoop.LatticePopulation(N)

transaction_frequencies = zeros(5) #the "new" way of allocating transactions
transaction_frequencies_old = zeros(5) # the "old" way of allocating transactions

for j in 1:numruns

    # track each set of neighbors, sorted
    doublets = []
    for (i, indv) in enumerate(CartesianIndices(pop.lattice))
        for (ni, neighbor) in enumerate(pop.neighbors[indv])
            doublet = sort([indv, neighbor])
            if doublet ∉ doublets
                push!(doublets, doublet)
            end
        end
    end

    transactions = Dict{CartesianIndex, Int64}()
    # randomly make one individual in each neighbor set the initiator
    # track the number of transactions each individual initiates
    for doublet in doublets
        initiator = doublet[rand(Bool)+1]
        if initiator ∉ keys(transactions)
            transactions[initiator] = 0
        end
        transactions[initiator] += 1
    end

    transactions_old = Dict{CartesianIndex, Int64}()
    ordering = LatticeCoop.assign_transactions(pop)
    # randomly order individuals

    # tmp_payoffs = LatticeCoop.get_payoffs(pop, ordering)
    # for indv in CartesianIndices(tmp_payoffs)
    #     transactions_2[indv] = -round(Int64, tmp_payoffs[indv])
    # end

    # the commented out block actually does the same thing:
    # check the "ordering" of each individual and assign transactions thereby
    for indv in CartesianIndices(ordering)
        for neighbor in pop.neighbors[indv]
            if indv ∉ keys(transactions_old)
                transactions_old[indv] = 0
            end
            transactions_old[indv] += ordering[indv] > ordering[neighbor]
        end
    end


    # hacked together way of counting the frequency of transaction initiations
    for indv in CartesianIndices(ordering)
        # N^2 because that's how many neighbor pairs there are
        transaction_frequencies_old[transactions_old[indv]+1] += 1.0/N^2/numruns
        if indv ∉ keys(transactions)
            transaction_frequencies[1] += 1.0/N^2/numruns
        else
            transaction_frequencies[transactions[indv]+1] += 1.0/N^2/numruns
        end
    end
end

fig = plt.figure()
plt.plot(transaction_frequencies, label="new method")
plt.plot(transaction_frequencies_old, label="old method")
ax=plt.gca()
ax.legend(loc=2)
ax.set_ylabel("frequency")
ax.set_xlabel("transactions performed by individual")
plt.tight_layout()
plt.gcf()

# sanity check: transaction frequencies should sum to 1
println("sum of \"new\" frequencies: $(sum(transaction_frequencies))")
println("sum of \"old\" frequencies: $(sum(transaction_frequencies_old))")
