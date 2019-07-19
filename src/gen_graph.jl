#!/usr/bin/env julia

## gen_graph.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test methods for generating random regular graphs.

using Random, StatsBase, Combinatorics, Revise
import LightGraphs.random_regular_graph
using LightGraphs:getRNG, SimpleEdge

N = 6
k = 3
g1 = 3
g2 = 0

rng = getRNG(-1)

function _suitable(edges::Set{SimpleEdge{T}}, potential_edges::Dict{T,T}) where T <: Integer
    println("checking suitability")
    println("edges = $edges, potential_edges = $potential_edges")
    isempty(potential_edges) && return true
    list = keys(potential_edges)
    for s1 in list, s2 in list
        s1 >= s2 && continue
        (SimpleEdge(s1, s2) ∉ edges) && return true
    end
    return false
end

_try_creation(n::Integer, k::Integer, rng::AbstractRNG) = _try_creation(n, fill(k, n), rng)

function _try_creation(n::T, k::Vector{T}, rng::AbstractRNG) where T <: Integer
    edges = Set{SimpleEdge{T}}()
    m = 0
    stubs = zeros(T, sum(k))
    for i = one(T):n
        for j = one(T):k[i]
            m += 1
            stubs[m] = i
        end
    end
    println("stubs = $stubs")
    # stubs = vcat([fill(i, k[i]) for i = 1:n]...) # slower

    while !isempty(stubs)
        potential_edges =  Dict{T,T}()
        println("potential_edges = $potential_edges")
        shuffle!(rng, stubs)
        for i in 1:2:length(stubs)
            s1, s2 = stubs[i:(i + 1)]
            println("$i, $s1, $s2")
            if (s1 > s2)
                s1, s2 = s2, s1
            end
            e = SimpleEdge(s1, s2)
            if s1 != s2 && ∉(e, edges)
                push!(edges, e)
            else
                potential_edges[s1] = get(potential_edges, s1, 0) + 1
                potential_edges[s2] = get(potential_edges, s2, 0) + 1
            end
            println("edges = $edges, potential_edges = $potential_edges")
        end

        if !_suitable(edges, potential_edges)
            return Set{SimpleEdge{T}}()
        end

        stubs = Vector{Int}()
        for (e, ct) in potential_edges
            append!(stubs, fill(e, ct))
        end
    end
    return edges
end

function random_regular_graph(n::Integer, k::Integer; seed::Int=-1)
    !iseven(n * k) && throw(ArgumentError("n * k must be even"))
    !(0 <= k < n) && throw(ArgumentError("the 0 <= k < n inequality must be satisfied"))
    if k == 0
        return SimpleGraph(n)
    end
    if (k > n / 2) && iseven(n * (n - k - 1))
        return complement(random_regular_graph(n, n - k - 1, seed=seed))
    end

    rng = getRNG(seed)

    edges = _try_creation(n, k, rng)
    #println("edges = $edges")
    while isempty(edges)
        edges = _try_creation(n, k, rng)
    end

    g = SimpleGraph(n)
    for edge in edges
        add_edge!(g, edge)
    end

    return g
end

gg = random_regular_graph(N,k)
println("$gg")
#
# g1_neighbors = [Int64[] for i in 1:N]
# if g1 == 2
#     for i in 1:N
#         j = 1
#         if i+j > N
#             push!(g1_neighbors[i], i+j-N)
#         else
#             push!(g1_neighbors[i], i+j)
#         end
#         if i-j < 1
#             push!(g1_neighbors[i], i-j+N)
#         else
#             push!(g1_neighbors[i], i-j)
#         end
#     end
# end
#
# unfinished_indvs = Int64[]
# if g1 > 2
#     for i in 1:N
#         for j in 1:N
#             println("$i, $j")
#             if length(g1_neighbors[j]) < g1
#                 push!(unfinished_indvs, j)
#             end
#         end
#         if length(unfinished_indvs) == 0
#             break
#         end
#         host = unfinished_indvs[1]
#         deleteat!(unfinished_indvs, 1)
#         while g1 - length(g1_neighbors[host]) > 0
#             if length(unfinished_indvs) == 0
#                 break
#             end
#             guest = rand(1:length(unfinished_indvs))
#             println("$i, $host, $guest, $unfinished_indvs")
#             push!(g1_neighbors[host], unfinished_indvs[guest])
#             push!(g1_neighbors[unfinished_indvs[guest]], host)
#             deleteat!(unfinished_indvs, guest)
#         end
#         unfinished_indvs == []
#     end
# end
#
# println("$g1_neighbors")

#
# function isregular(
#     neighbor_list::Array{Array{Int64,1}},
#     degree::Int64
#     )
#     num_appearances = zeros(Int64, length(neighbor_list))
#     #println("$neighbor_list")
#     for (li, list) in enumerate(neighbor_list)
#         if length(list) != degree
#             return false
#         elseif li ∈ list
#             return false
#         end
#         for ii in list
#             num_appearances[ii] += 1
#         end
#     end
#     #println("$num_appearances")
#     if any(num_appearances .!= degree)
#         return false
#     else
#         return true
#     end
# end
#
# num_attempts = 0
# successes = 0
#
# g1_neighbors = [Int64[] for x in 1:N]
# g2_neighbors = [Int64[] for x in 1:N]
# all_neighbors = [Int64[] for x in 1:N]
#
# S = collect(combinations(1:N,g1))
#
# while !isregular(g1_neighbors, g1)
#     global num_attempts += 1
#     [g1_neighbors[x] = Int64[] for x in 1:N]
#     S = collect(combinations(1:N,2))
#     while !isempty(S)
#         #println(S)
#         wts = zeros(length(S))
#         for (pi, (i,j)) in enumerate(S)
#             wts[pi] = (g1 - length(g1_neighbors[i]))*
#                 (g1 - length(g1_neighbors[j]))
#         end
#         #newpair = S[sample(1:length(S), Weights(wts))]
#         newpair = sample(S, Weights(wts))
#         push!(g1_neighbors[newpair[1]], newpair[2])
#         push!(g1_neighbors[newpair[2]], newpair[1])
#         filter!(x -> x != newpair, S)
#         for (pi, (i,j)) in enumerate(S)
#             if length(g1_neighbors[i]) > (g1 - 1) ||
#                 length(g1_neighbors[j]) > (g1 - 1)
#                 filter!(x -> x != [i,j], S)
#             end
#         end
#     end
#     #println("$g1_neighbors")
# end
# println("$g1_neighbors")
# g1_edges = Array{Int64,1}[]
# for i in 1:N
#     for j in 1:g1
#         edge = [i, g1_neighbors[i][j]]
#         if issorted(edge)
#             push!(g1_edges, edge)
#         end
#     end
# end
# all_neighbors = g1_neighbors
#
# while !isregular(g2_neighbors, g2)
#     global num_attempts += 1
#     [g2_neighbors[x] = Int64[] for x in 1:N]
#     S = collect(combinations(1:N,2))
#     filter!(x -> !(x ∈ g1_edges), S)
#     while !isempty(S)
#         #println(S)
#         wts = zeros(length(S))
#         for (pi, (i,j)) in enumerate(S)
#             #println("$pi, $i, $j")
#             #println("$(g1 - length(g1_neighbors[i])), $(g1 - length(g1_neighbors[j]))")
#             wts[pi] = (g2 - length(g2_neighbors[i]))*
#                 (g2 - length(g2_neighbors[j]))
#         end
#         #newpair = S[sample(1:length(S), Weights(wts))]
#         newpair = sample(S, Weights(wts))
#         push!(g2_neighbors[newpair[1]], newpair[2])
#         push!(g2_neighbors[newpair[2]], newpair[1])
#         filter!(x -> x != newpair, S)
#         for (pi, (i,j)) in enumerate(S)
#             if length(g2_neighbors[i]) > (g2 - 1) ||
#                 length(g2_neighbors[j]) > (g2 - 1)
#                 filter!(x -> x != [i,j], S)
#             end
#         end
#     end
#     #println("$g2_neighbors")
# end
# println("$g2_neighbors")
# all_neighbors = [vcat(all_neighbors[x], g2_neighbors[x]) for x in 1:N]
# println("$num_attempts")

# for Ni in 1:(N-1)
#     unfilled_indvs = filter(x -> !(x ∈ filled_indvs_g1),
#         (Ni+1):N)
#     indvs_to_add = sample(unfilled_indvs,
#         g1-length(g1_neighbors[Ni]), replace=false)
#     println("$Ni, $unfilled_indvs, $indvs_to_add, $g1_neighbors")
#     [push!(g1_neighbors[Ni], ii) for ii in indvs_to_add]
#     [push!(g1_neighbors[ii], Ni) for ii in indvs_to_add]
#     [push!(all_neighbors[Ni], ii) for ii in indvs_to_add]
#     [push!(all_neighbors[ii], Ni) for ii in indvs_to_add]
#     push!(filled_indvs_g1, Ni)
#     for ii in indvs_to_add
#         if length(g1_neighbors[ii]) == g1
#             push!(filled_indvs_g1, ii)
#         end
#     end
# end
