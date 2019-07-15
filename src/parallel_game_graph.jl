#!/usr/bin/env julia

## test_game_graph.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Implement and test GraphGame.

using Random, LightGraphs
using Distributed
using Revise
using GraphGame
using ArgParse
using CSV
using DataFrames
import JSON



#rc("font", size="large")

function read_parameters(
    defpars::Dict{String, Any},
    inputfile=nothing
    )

    pars=copy(defpars)

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
            convertf = (val) -> round(T, val)
        else
            convertf = (val) -> convert(T, val)
        end

        # use defpars for list of usable parameters in JSON
        if parkey in keys(inpars)
            if "value" in keys(inpars[parkey])
                val = inpars[parkey]["value"]
            elseif "range" in keys(inpars[parkey])
                valr = inpars[parkey]["range"]
                if "log" in keys(valr)
                    b = valr["log"]
                    rf = (r) -> b .^ r
                    pop!(valr, "log")
                else
                    rf = (r) -> r
                end
                start = pop!(valr, "start")
                rkws = Dict(zip(Symbol.(keys(valr)), values(valr)))
                val = rf(range(start; rkws...))
            end
        else
            val = pars[parkey]["value"]
        end

        if !isstructtype(typeof(val)) || typeof(val) == String
            pars[parkey] = [convertf(val)]
        else
            pars[parkey] = convertf.(val)
        end
    end

    return pars
end

function main(
    args
    )

    s = ArgParseSettings(description =
        "run GameGraphs simulation across multiple cores")
    @add_arg_table s begin
        "--ncpus"
            arg_type=Int64
            default=max(round(Int,Sys.CPU_THREADS/2), 1)
        "--input"
            default=nothing
        #"--output"
        #    default=nothing
    end
    parsed_args = parse_args(args, s)

    defpars = Dict{String,Any}([
        "N"     => Dict("value" => 10000, "type" => Int64),
        "k"     => Dict("value" => 3,     "type" => Int64),
        "w"     => Dict("value" => 1e-3,  "type" => Float64),
        "R"     => Dict("value" => 5,     "type" => Int64),
        "S"     => Dict("value" => 0,     "type" => Int64),
        "T"     => Dict("value" => 8,     "type" => Int64),
        "P"     => Dict("value" => 1,     "type" => Int64),
        "init_freq" => Dict("value" => 0.5, "type" => Float64),
        "update_type" => Dict("value" => "bd", "type" => String),
        "num_gens" => Dict("value" => 100, "type" => Int64),
        "num_trials" => Dict("value" => 1, "type" => Int64),
        "output" => Dict("value" => "test.csv", "type" => String)
    ])
    pars = read_parameters(defpars, parsed_args["input"])

    # take the Cartesian product of all parameter combinations
    parsets = collect(Base.product(values(pars)...))
    nsets = length(parsets)

    # setup workers assuming directory is manually added to LOAD_PATH
    addprocs(min(parsed_args["ncpus"], round(Int64, Sys.CPU_THREADS/2)))
    wpool = WorkerPool(workers())
    extradir = filter((p)->match(r"/", p) !== nothing, LOAD_PATH)[1]
    @everywhere workers() push!(LOAD_PATH, $extradir)
    @everywhere workers() eval(:(using Random))
    @everywhere workers() eval(:(using LightGraphs))
    @everywhere workers() eval(:(using GraphGame))
    @everywhere workers() eval(:(using Dates))

    inputs  = RemoteChannel(()->Channel{Dict}(2*nsets*maximum(pars["num_trials"])))
    results = RemoteChannel(()->Channel{Dict}(2*nsets*maximum(pars["num_trials"])))

    @everywhere function run_worker(inputs, results)
        # save trial number and random seed
        seed = Dict(zip(["seed1", "seed2", "seed3", "seed4"], Random.GLOBAL_RNG.seed))

        while true
            pard = take!(inputs)
            pard = merge(pard, seed)

            R, S, T, P = pard["R"], pard["S"], pard["T"], pard["P"]
            N, degree, w = pard["N"], pard["k"], pard["w"]
            update_type = pard["update_type"]
            init_freq = pard["init_freq"]
            num_gens, num_trials = pard["num_gens"], pard["num_trials"]

            println("--- running ", pard["nrun"], " --- ")
            flush(stdout)

            game = Float64[R S; T P]
            graph = random_regular_graph(N, degree)

            start = now()

            gg = GameGraph(N,
                degree,
                w,
                game,
                graph
                )
            num_defectors = floor(Int64, N*(1-init_freq))
            indvs_to_flip = randperm(N)[1:num_defectors]
            gg.strategies[indvs_to_flip] = ones(Int64, num_defectors)*2
            initialize_fitnesses!(gg)

            freq_trajectory = Array{Float64, 1}[]
            pair_freq_trajectory = Array{Float64, 1}[]
            cond_freq_trajectory = Array{Float64, 1}[]

            for i in 1:num_gens
                #println(i)
                freq_stats = get_conditional_frequencies(gg, true)
                push!(freq_trajectory, freq_stats[1])
                push!(pair_freq_trajectory, vec(freq_stats[2]))
                push!(cond_freq_trajectory, vec(freq_stats[3]))
                evolve!(gg, N, update_type)
                #println("$(gg.generation), $(get_frequencies(gg)), $(get_pair_frequencies(gg))")
            end

            # save frequency trajectories
            freq_trajectory = hcat(freq_trajectory...)
            pair_freq_trajectory = hcat(pair_freq_trajectory...)
            pard["p_c_trajectory"] = freq_trajectory[1,:]
            pard["p_cc_trajectory"] = pair_freq_trajectory[1,:]

            # output elapsed time
            stop = now()
            println("--- ran ", pard["nrun"], " --- elapsed time: ",
                Dates.canonicalize(Dates.CompoundPeriod(round(stop-start, Dates.Second(1)))))
            foreach(k -> print(k, ": ", pard[k], ", "),
                sort(collect(keys(filter(p->p.first ∉ ["p_c_trajectory", "p_cc_trajectory"], pard)))))
            println()
            flush(stdout)

            # return data to master process
            put!(results, pard)
        end
    end

    # load parameter sets into inputs channel
    nruns = 0
    for parset in parsets
        pard = Dict(zip(keys(pars), parset))
        #println(pard)
        println("--- queueing --- ")
        foreach(k -> print(k, ": ", pard[k], ", "), sort(collect(keys(pard))))
        println()
        flush(stdout)
        for rep in 1:pard["num_trials"]
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
    file = occursin(r"\.csv$", output) ? output : output*".csv"
    cols = push!(sort(collect(keys(pars))),
                 ["rep", "p_c_trajectory", "p_cc_trajectory", "seed1", "seed2", "seed3", "seed4"]...)
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
end

#main(ARGS)

main(["--input", "submit/ohtsuki_PD_1.json"])
