#!/usr/bin/env julia

## test_game_graph.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Implement and test GraphGame.

using Random, LightGraphs
using Distributed
using Revise
using ArgParse
using CSV
using Dates
using DataFrames
import JSON



#rc("font", size="large")

function read_parameters(defpars::Dict{String,Any},
    inputfile = nothing)

    pars = copy(defpars)

    # read JSON file
    if inputfile != nothing
        inpars = JSON.parsefile(inputfile)
    else
        inpars = Dict()
    end

    for parkey in keys(defpars)
        if "type" in keys(pars[parkey])
            if isprimitivetype(pars[parkey]["type"]) ||
                pars[parkey]["type"] == String
                T = pars[parkey]["type"]
            end
        else
            # default type is Float64
            T = Float64
        end
        #println(parkey, T)
        if T <: Int
            convertf = (val)->round(T, val)
        else
            convertf = (val)->convert(T, val)
        end

        # use defpars for list of usable parameters in JSON
        if parkey in keys(inpars)
            if "value" in keys(inpars[parkey])
                val = inpars[parkey]["value"]
            elseif "range" in keys(inpars[parkey])
                valr = inpars[parkey]["range"]
                if "log" in keys(valr)
                    b = valr["log"]
                    rf = (r)->b.^r
                    pop!(valr, "log")
                else
                    rf = (r)->r
                end
                start = pop!(valr, "start")
                rkws = Dict(zip(Symbol.(keys(valr)), values(valr)))
                val = rf(range(start; rkws...))
            end
        else
            val = pars[parkey]["value"]
        end

        if !isstructtype(typeof(val)) || typeof(val) == String || typeof(val) == Bool
            pars[parkey] = [convertf(val)]
        else
            pars[parkey] = convertf.(val)
        end
    end

    return pars
end

function main(args)

    s = ArgParseSettings(description =
        "run GameGraphs simulation across multiple cores")
    @add_arg_table s begin
        "--ncpus"
            arg_type = Int64
            default = max(round(Int, Sys.CPU_THREADS / 2), 1)
        "--input"
            default = nothing
        #"--output"
        #    default=nothing
    end
    parsed_args = parse_args(args, s)

    defpars = Dict{String,Any}([
        "N"     => Dict("value" => 200, "type" => Int64),
        "k"     => Dict("value" => 6, "type" => Int64),
        "g1"     => Dict("value" => 1,     "type" => Int64),
        "test_strat" => Dict("value" => 1, "type" => Int64),
        "w"     => Dict("value" => 1e-2,  "type" => Float64),
        "BC_ratio" => Dict("value" => 1.0, "type" => Float64),
        "update_type" => Dict("value" => "db", "type" => String),
        "fitness_function" => Dict("value" => "DoL_add", "type" => String),
        "threshold" => Dict("value" => 1, "type" => Int64),
        "num_trials" => Dict("value" => 100, "type" => Int64),
        "num_runs" => Dict("value" => 100, "type" => Int64),
        "runs_per_graph" => Dict("value" => 10, "type" => Int64),
        "output" => Dict("value" => "output/test.csv", "type" => String)
    ])
    pars = read_parameters(defpars, parsed_args["input"])

    # take the Cartesian product of all parameter combinations
    parsets = collect(Base.product(values(pars)...))
    nsets = length(parsets)

    # setup workers assuming directory is manually added to LOAD_PATH
    addprocs(min(parsed_args["ncpus"], round(Int64, Sys.CPU_THREADS / 2)))
    wpool = WorkerPool(workers())
    extradir = filter((p)->match(r"/", p) !== nothing, LOAD_PATH)[1]
    @everywhere workers() push!(LOAD_PATH, $extradir)
    @everywhere workers() eval(:(using Random))
    @everywhere workers() eval(:(using GenGraph))
    @everywhere workers() eval(:(using MultiEdgeGame))
    @everywhere workers() eval(:(using Dates))

    inputs  = RemoteChannel(()->Channel{Dict}(2 * nsets * maximum(pars["num_trials"])))
    results = RemoteChannel(()->Channel{Dict}(2 * nsets * maximum(pars["num_trials"])))

    @everywhere function run_worker(inputs, results)
        # save trial number and random seed
        seed = Dict(zip(["seed1", "seed2", "seed3", "seed4"], Random.GLOBAL_RNG.seed))

        while true
            pard = take!(inputs)
            pard = merge(pard, seed)

            BC_ratio = pard["BC_ratio"]
            w = pard["w"]
            n, g1, g2 = pard["N"], pard["g1"], pard["k"] - pard["g1"]
            update_type = pard["update_type"]
            num_runs, num_trials = pard["num_runs"], pard["num_trials"]
            runs_per_graph = pard["runs_per_graph"]
            fitness_function = pard["fitness_function"]
            threshold = pard["threshold"]
            test_strat = pard["test_strat"]

            println("--- running ", pard["nrun"], " --- ")
            flush(stdout)

            game = get_game_function(fitness_function)

            g = [g1, g2]

            C = 1
            B = BC_ratio

            fixed = 0
            extinct = 0

            start = now()

            for i in 1:num_trials
                graph = generate_multitype_connected_graph(n, g)
                for j in 1:runs_per_graph
                    pop = Population(n, ones(Int64, n) * (2 - (test_strat == 2)), zeros(n), 0)
                    mgg = MultiGameGraph(n, g, 2, graph[1], graph[2], game, w, GameParams(B, C, threshold))
                    pop.strategies[rand(1:n)] = test_strat
                    initialize_fitnesses!(pop, mgg)

                    freq = get_frequency(pop, mgg, test_strat)
                    while freq ∉ (0, 1)
                        evolve!(pop, mgg, n, update_type)
                        freq = get_frequency(pop, mgg, test_strat)
                    end
                    if freq == 1
                        fixed += 1
                    else
                        extinct += 1
                    end
                end
            end
            pard["fixed"] = fixed
            pard["extinct"] = extinct
            stop = now()
            println("--- ran ", pard["nrun"], " --- elapsed time: ",
                Dates.canonicalize(Dates.CompoundPeriod(round(stop - start, Dates.Second(1)))))
            foreach(k->print(k, ": ", pard[k], ", "),
                sort(collect(keys(filter(p->p.first ∉ ["nrun"], pard)))))
            println()
            flush(stdout)

            # return data to master process
            put!(results, pard)
        end
    end

    total_time_start = now()

    # load parameter sets into inputs channel
    nruns = 0
    for parset in parsets
        pard = Dict(zip(keys(pars), parset))
        #println(pard)
        println("--- queueing --- ")
        foreach(k->print(k, ": ", pard[k], ", "), sort(collect(keys(pard))))
        println()
        flush(stdout)
        for rep in 1:pard["num_runs"]
            nruns += 1
            rpard = copy(pard)
            rpard["rep"] = rep
            rpard["nrun"] = nruns
            put!(inputs, rpard)
        end
    end

    # start workers running on parameter sets in inputs
    for w in workers() # start tasks on the workers to process requests in parallel
        remote_do(run_worker, w, inputs, results)
    end

    # create output file name and data table
    output = pars["output"][1]
    println(output)
    file = occursin(r"\.csv$", output) ? output : output * ".csv"
    cols = push!(sort(collect(keys(pars))),
                 ["rep", "fixed", "extinct", "seed1", "seed2", "seed3", "seed4"]...)
    dat = DataFrame(Dict([(c, Any[]) for c in cols]))

    # grab results and output to CSV
    for sim in 1:nruns
        # get results from parallel jobs
        flush(stdout)
        resd = take!(results)
        nrun = pop!(resd, "nrun")

        # add to table (must convert dict keys to symbols) and save
        push!(dat, Dict([(Symbol(k), resd[k]) for k in keys(resd)]))
        CSV.write(file, dat)
    end
    total_time_stop = now()
    total_time = Dates.canonicalize(Dates.CompoundPeriod(round(total_time_stop - total_time_start, Dates.Second(1))))
    println("total time elapsed: $total_time")
end

#main(ARGS)

main(["--input", "submit/su_fig_threshold_2.json"])
