#!/usr/bin/env julia

## gen_graph.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test methods for generating random regular graphs.

# TODO: check and see if the graph is connected.
# If not, we will run into problems when trying to calculate
# fixation probabilities.

using Random, StatsBase, Combinatorics, Revise

N = 500
g = [1, 1]


function is_suitable(
    edges::Array{Array{Int64,1}},
    potential_edges::Dict{Int64,Int64},
    extant_edges::Array{Array{Int64,1}} = [Int64[]]
    )
    # given a Dict of potential edges,
    # see if there's any way to construct a new edge

    # if the list is empty, everything is hunky dory
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
    all_neighbors::Array{Array{Int64,1},1},
    curr_connections::Array{Int64,1},
    indv::Int64,
    )
    if indv ∉ curr_connections
        push!(curr_connections, indv)
    end
    for neighb in all_neighbors[indv]
        if neighb ∉ curr_connections
            push!(curr_connections, neighb)
            #println("$indv, $neighb, $(vcat(tmp_connections, curr_connections ...))")
            new_connections = assemble_connected_list(all_neighbors, curr_connections, neighb)
            #println("$new_connections")
            #new_connections = assemble_connected_list(all_neighbors, curr_connections, neighb)
            [push!(curr_connections, x) for x in new_connections if x ∉ curr_connections]
        end
    end
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
    let
    connected = false
    while !(connected)
        global gg = generate_multitype_random_regular(N, g)
        connected = is_connected(gg[2])
        num_attempts += 1
        println("$num_attempts")#, $(gg[2])")
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
