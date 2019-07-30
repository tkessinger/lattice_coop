#!/usr/bin/env julia

## MultiGraphGame.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test methods for generating random regular graphs.

# TODO: check and see if the graph is connected.
# If not, we will run into problems when trying to calculate
# fixation probabilities.

using Random, StatsBase, Combinatorics, Revise

N = 500
g = [1, 1]

mutable struct MultiGameGraph
    # by design this is a random regular graph

    n::Int64 # number of nodes
    g::Array{Int64, 1} # degrees of edge-typed subgraphs
    w::Float64 # strength of selection
    game::Function # an associated game to be played on the graph
    graph::Dict{Int64, Array{Array{Int64, 1}, 1}} # a Dict of individuals and their neighbors by type

    # the next three are the only mutable attributes
    generation::Int64 # current generation
    strategies::Array{Int64, 1} # an array of individuals' strategies
    fitnesses::Array{Float64, 1} # individuals' fitnesses depending on,
        # e.g., neighbor payoffs
end

function DoL_payoff(
    graph::Dict{Int64, Array{Int64, 2}},
    strategies::Array{Int64, 1},
    indv::Int64
    )
    s1 = sum([strategies[x] .== 1 for x in graph[indv][1]]) # number of cooperators along edgetype 1
    s2 = sum([strategies[x] .== 1 for x in graph[indv][2]]) # number of cooperators along edgetype 2
    # if indv is a cooperator
    if strategies[indv] == 1
        if s2 >= 1
            return (s1 + s2 + 1)*(B-C)
        else
            return -C
        end
    # if indv is a defector
    elseif strategies[indv] == 2
        if s1 >= 1 && s2 >= 1
            return (s1 + s2)*B
        else
            return 0
        end
    end
end

function initialize_fitnesses!(
    gg::MultiGameGraph
    )
    # assign the correct fitness values to each individual in the population
    # ideally this would be done within a constructor
    for indv in keys(gg.graph)
        payoff = gg.game(gg.graph, gg.strategies, indv)
        fitness = 1 - gg.w + gg.w * payoff
        gg.fitnesses[indv] = fitness
    end
end

function update_indv!(
    gg::MultiGameGraph,
    invader::Int64,
    invadee::Int64
    )
    # if an individual gets replaced, we need to update their fitness
    # as well as that of all their neighbors
    update_strategies!(gg, invader, invadee)
    update_fitness!(gg, invadee)
end

function update_strategies!(
    gg::MultiGameGraph,
    invader::Int64,
    invadee::Int64
    )
    # replaces the strategy of a replaced individual
    gg.strategies[invadee] = gg.strategies[invader]
end

function update_fitness!(
    gg::MultiGameGraph,
    indv::Int64
    )
    # updates the fitness of an individual whose strategy is replaced
    # as well as the fitnesses of all their neighbors
    payoff = gg.game(gg.graph, gg.strategies, indv)
    fitness = 1.0 - gg.w + gg.w * payoff
    return fitness
end

function evolve!(
    gg::MultiGameGraph,
    num_gens::Int64=1,
    update_rule::String="birth_death"
    )
    # depending on the update rule and their fitnesses,
    # replaces an individual with a neighbor

    for i in 1:num_gens
        if update_rule == "birth_death" || update_rule == "bd"
            # an individual is chosen to reproduce proportional to fitness
            # they replace a random neighbor
            indv = sample(1:gg.n, Weights(gg.fitnesses))
            neighbor_to_replace = rand(neighbors(gg.graph, indv))
            # this if block is added in case there's no need to update
            if gg.strategies[indv] != gg.strategies[neighbor_to_replace]
                update_indv!(gg, indv, neighbor_to_replace)
            end
        elseif update_rule == "death_birth" || update_rule == "db"
            # an individual is chosen to die at random
            # their neighbors compete to replace them, proportional to fitness
            indv = rand(1:gg.n)
            neighbor_fitnesses = gg.fitnesses[neighbors(gg.graph, indv)]
            invading_neighbor = sample(neighbors(gg.graph, indv), Weights(neighbor_fitnesses))
            if gg.strategies[invading_neighbor] != gg.strategies[indv]
                update_indv!(gg, invading_neighbor, indv)
            end
        elseif update_rule == "imitation" || update_rule == "im"
            # an individual compares its fitness against its neighbors and itself
            # then adopts their strategy, proportional to fitness
            indv = rand(1:gg.n)
            neighbor_fitnesses = gg.fitnesses[vcat(indv, neighbors(gg.graph, indv))]
            invading_neighbor = sample(vcat(indv, neighbors(gg.graph, indv)), Weights(neighbor_fitnesses))
            if gg.strategies[invading_neighbor] != gg.strategies[indv]
                update_indv!(gg, invading_neighbor, indv)
            end
        end
    gg.generation += 1
    end
end

function is_suitable(
    edges::Array{Array{Int64,1}},
    potential_edges::Dict{Int64,Int64},
    extant_edges::Array{Array{Int64,1}} = [Int64[]]
    )
    # given a Dict of potential edges,
    # see if there's any way to construct a new edge

    # if the Dict is empty, everything is hunky dory
    if isempty(potential_edges)
        return true
    end
    list = keys(potential_edges)
    # if it's not empty, check if a new edge can be created
    for s1 in list, s2 in list
        s1 >= s2 && continue
        # if an edge doesn't already exist, that's good
        if [s1, s2] ∉ edges && [s1, s2] ∉ extant_edges
            return true
        end
    end
    # if there's no way to create a new edge, return false
    return false
end

function try_random_regular(
    N::Int64,
    k::Int64,
    extant_edges::Array{Array{Int64,1},1} = Array{Int64,1}[]
    )
    # attempt to generate a random regular graph
    # based on a degree and a list of current edges
    # if this fails, return an empty list: the ``main'' function will try again


    # list of ``unfilled'' edges
    stubs = vcat([fill(i,k) for i in 1:N] ...)
    # the empty edge list of the current type, which we will populate
    edges = Array{Int64,1}[]
    # if the stub list isn't empty, let's try to empty it
    while !isempty(stubs)
        # a Dict() containing the number of possible edges
        # that each stub could be part of
        potential_edges = Dict{Int64,Int64}()
        shuffle!(stubs)
        # consider random pairs of stubs
        for i in 1:2:length(stubs)
            s1, s2 = stubs[i:(i+1)]
            # order the stubs
            if s1 > s2
                s1, s2 = s2, s1
            end
            # if the edge isn't a loop and doesn't already appear anywhere,
            # then add it to the list of edges
            if s1 != s2 && ∉([s1, s2], edges) && ∉([s1, s2], extant_edges)
                # println("adding edge [$s1, $s2]")
                push!(edges, [s1, s2])
            # if this condition fails, add 1 to potential_edges
            else
                potential_edges[s1] = get(potential_edges, s1, 0) + 1
                potential_edges[s2] = get(potential_edges, s2, 0) + 1
            end
        end
        # see if there's any way to construct a new edge based on what we've got
        if !(is_suitable(edges, potential_edges, extant_edges))
            # println("not suitable!")
            # if not, return an empty edge list
            return Array{Int64, 1}[]
        end

        # repopulate the stub list based on what remains, then try again
        stubs = Int64[]
        for (e, ct) in potential_edges
            append!(stubs, fill(e, ct))
        end
        #println("$stubs, $edges, $potential_edges")
    end
    return edges
end

function generate_random_regular(
    N::Int64,
    k::Int64
    )

    num_attempts = 0

    neighbors = [Int64[] for i in 1:N]

    edges = []
    while isempty(edges)
        num_attempts += 1
        edges = try_random_regular(N, k)
    end
    for (s1, s2) in edges
        push!(neighbors[s1], s2)
        push!(neighbors[s2], s1)
    end
    println("$num_attempts")
    return neighbors
end

function generate_multitype_random_regular(
    N::Int64,
    k::Int64
    )
    return generate_multitype_random_regular(N, [k])
end

function generate_multitype_random_regular(
    N::Int64, # number of individuals
    g::Array{Int64, 1} # list of number of edges of each type
    )
    # generates a graph (list of neighbors for each individual)
    # in which there are multiple edge types

    # sum(g) needs to be less than N
    if sum(g) > N
        throw(DomainError(g, "sum(g) must be <= N"))
    end

    # list of edges by type
    # edges[i][j] should give the jth edge of type i
    edges = [Array{Int64,1}[] for x in 1:length(g)]
    # println("$edges")
    # list of neighbors by individual and types
    # neighbors[i][j] should give the neighbors of type i for individual j
    neighbors = [[Int64[] for i in 1:N] for x in 1:length(g)]

    # if any type of edge list is empty, try re-generating the graph
    # this function will return empty edge lists if there's no way to make them work
    while any([isempty(edgelist) for edgelist in edges])
        edges = recursive_graph_gen(N, g, length(g))
    end
    # populate neighbors based on edges
    for (i, gi) in enumerate(g)
        for (s1, s2) in edges[i]
            push!(neighbors[i][s1], s2)
            push!(neighbors[i][s2], s1)
        end
    end

    all_neighbors = [Int64[] for i in 1:N]
    [[[push!(all_neighbors[indv], neighbor) for neighbor in neighbors[ii][indv]] for indv in 1:N] for ii in 1:length(g)]

    # sort the neighbor list to make it easier to read
    [sort!(neighbs) for neighbs in all_neighbors]
    [[sort!(neighbs) for neighbs in neighbors[i]] for i in 1:length(g)]
    return neighbors, all_neighbors
end

function recursive_graph_gen(
    N::Int64, # number of individuals
    g::Array{Int64, 1}, # number of edges of each type
    ii::Int64 # current type of edge we are trying to generate
    )
    # attempt to generate a multitype-edge regular graph
    # if this turns out to be impossible
    # (e.g., because one type of edge is assigned in a way that
    # makes other types impossible)
    # then return an empty list for that type

    edges = Array{Array{Int64,1},1}[]
    # recursively call this function if we haven't ``hit'' the lowest level yet
    if ii > 1
        edges = recursive_graph_gen(N, g, (ii-1))
    end
    # ugly hack needed to make extant_edges have the right type
    extant_edges = Array{Array{Int64,1},1}(vcat(edges ...))
    # generate a list of edges based on the ``lower'' levels
    tmp_edges = try_random_regular(N, g[ii], extant_edges)
    #println("$edges, $tmp_edges, $ii")
    push!(edges, tmp_edges)
    return edges
end

function assemble_connected_list(
    all_neighbors::Array{Array{Int64,1},1}, # list of neighbors of each individual
    curr_connections::Array{Int64,1}, # list containing individuals that have already been counted
    indv::Int64, # current individual
    )
    # create a list of all individuals the current individual is connected to
    # this is used to determine whether the graph is connected or not

    # if the current individual isn't on the list, put them there
    if indv ∉ curr_connections
        push!(curr_connections, indv)
    end
    for neighb in all_neighbors[indv]
        if neighb ∉ curr_connections
            # add the neighbor to the list
            push!(curr_connections, neighb)
            #println("$indv, $neighb, $(vcat(tmp_connections, curr_connections ...))")
            # add all their connections, and all their connections (recursively)
            new_connections = assemble_connected_list(all_neighbors, curr_connections, neighb)
            #println("$new_connections")
            #new_connections = assemble_connected_list(all_neighbors, curr_connections, neighb)
            [push!(curr_connections, x) for x in new_connections if x ∉ curr_connections]
        end
    end
    # concatenate the list so ∈ will work
    return vcat(curr_connections ...)
end

function is_connected(
    neighbors::Array{Array{Array{Int64,1},1},1}
    )
    all_neighbors = [Int64[] for i in 1:N]
    [[[push!(all_neighbors[indv], neighbor) for neighbor in neighbors[ii][indv]] for indv in 1:N] for ii in 1:length(neighbors)]
    return is_connected(all_neighbors)
end

function is_connected(
    all_neighbors::Array{Array{Int64,1},1},
    )
    N = length(all_neighbors)
    connected_list = assemble_connected_list(all_neighbors, Int64[], 1)
    connected_list = sort(unique(connected_list))
    #println("$connected_list")
    if length(connected_list) == N
    #if connected_list == collect(1:N)
        return true
    else
        return false
    end
end

function generate_multitype_connected_graph(
    N::Int64, # number of individuals
    g::Array{Int64, 1} # list of number of edges of each type
    )
    num_attempts = 1
    let connected = false
        while !(connected)
            global gg = generate_multitype_random_regular(N, g)
            connected = is_connected(gg[2])
            num_attempts += 1
            #println("$num_attempts")#, $(gg[2])")
        end
    end
    return gg
end

gg = generate_multitype_connected_graph(N, g)
#gg = generate_multitype_random_regular(N, g)
#is_connected(gg[2])
#assemble_connected_list([[2,3],[1,3],[1,2]], Int64[], 1)
#assemble_connected_list(gg[2], Int64[], 1)

#println("$edges")
