#!/usr/bin/env julia

## gen_graph.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test methods for generating random regular graphs.

using Random, StatsBase, Combinatorics, Revise

N = 6
k = 4
g1 = 2
g2 = 2


function is_suitable(
    edges::Array{Array{Int64,1}},
    potential_edges::Dict{Int64,Int64},
    extant_edges::Array{Array{Int64,1}} = [Int64[]]
    )
    if isempty(potential_edges)
        return true
    end
    list = keys(potential_edges)
    for s1 in list, s2 in list
        s1 >= s2 && continue
        if [s1, s2] ∉ edges
            return true
        end
    end
    # if all([sort([s1,s2]) ∉ edges for s1 in list, s2 in list])
    #     return true
    # end
    return false
end

function try_random_regular(
    N::Int64,
    k::Int64,
    extant_edges::Array{Array{Int64,1}} = Array{Int64,1}[]
    )
    stubs = vcat([fill(i,k) for i in 1:N] ...)
    edges = Array{Int64,1}[]
    while !isempty(stubs)
        potential_edges = Dict{Int64,Int64}()
        shuffle!(stubs)
        for i in 1:2:length(stubs)
            s1, s2 = stubs[i:(i+1)]
            if s1 > s2
                s1, s2 = s2, s1
            end
            if s1 != s2 && ∉([s1, s2], edges) && ∉([s1, s2], extant_edges)
                push!(edges, [s1, s2])
            else
                potential_edges[s1] = get(potential_edges, s1, 0) + 1
                potential_edges[s2] = get(potential_edges, s2, 0) + 1
            end
        end
        println("$edges, $potential_edges")
        if !(is_suitable(edges, potential_edges))
            #println("not suitable!")
            return Array{Int64, 1}[]
        end

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
    g::Array{Int64, 1}
    )

    all_edges = Array{Int64,1}[]
    edges = [Array{Int64,1}[] for x in 1:length(g)]
    neighbors = [[Int64[] for i in 1:N] for x in 1:length(g)]

    for (i, gi) in enumerate(g)
        tmp_edges = []
        while isempty(tmp_edges)
            println("trying")
            tmp_edges = try_random_regular(N, gi, edges[1:(i-1)])
        end
        println("$tmp_edges")
        #[push!(all_edges, x) for x in tmp_edges]
        edges[i] = tmp_edges

    end
    for (i, gi) in enumerate(g)
        for (s1, s2) in edges[i]
            push!(neighbors[i][s1], s2)
            push!(neighbors[i][s2], s1)
        end
    return neighbors

end

gg = generate_multitype_random_regular(N, [g1, g2])

#println("$edges")
