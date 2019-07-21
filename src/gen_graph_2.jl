#!/usr/bin/env julia

## gen_graph.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test methods for generating random regular graphs.

# TODO: check and see if the graph is connected.
# If not, we will run into problems when trying to calculate
# fixation probabilities.

using Random, StatsBase, Combinatorics, Revise

N = 8
k = 6
g1 = 2
g2 = 2
g3 = 2


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

    # list of all edges (undirected)
    all_edges = Array{Int64,1}[]

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
    # sort the neighbor list to make it easier to read
    [[sort!(neighbs) for neighbs in neighbors[i]] for i in 1:length(g)]
    return neighbors

end

function recursive_graph_gen(
    N::Int64, # number of individuals
    g::Array{Int64, 1}, # number of edges of each type
    ii::Int64 # current type of edge we are trying to generate
    )
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

g = [g1, g2]
gg = generate_multitype_random_regular(N, [g1, g2, g3]);

#println("$edges")
