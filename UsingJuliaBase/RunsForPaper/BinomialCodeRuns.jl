import QuantumToolbox as qt
using Logging
using LoggingExtras

using Dates

using MiniLoggers

import SuperconductingCavities as SC

# MiniLogger(minlevel = MiniLoggers.Info) |> global_logger

batch_name = ARGS[2]
the_path = ARGS[3]
run_name = ARGS[4]

log_file = the_path*"/logs/"*run_name*".log"
#log_file2 = open(the_path*"/logs/"*run_name*"_progbars.log", "a")

# logger = FormatLogger(open(log_file, "a")) do io, args
#     df = DateFormat("e-u-d-yy:H:M")
#     t = now()
#     the_time = string(Dates.format(t, df))
#     println(io, the_time, ":", "[", args.level, "] ", args.message )
# end

# #logger = FileLogger(log_file; append = true)
# logger = MinLevelLogger(logger, LogLevel(-1))
# global_logger(logger)


log_file = open(log_file, "a")
InfoLogger = MiniLogger(io = log_file, minlevel = MiniLoggers.Info)
#ProgressLogger = MiniLogger(minlevel = LogLevel(-1))
#DebugLogger = MiniLogger(minlevel = MiniLoggers.Debug)

global_logger(InfoLogger)


models = Dict{Any, Any}("Mode3" => "ModelSaves/Mode3_ForBinomialBatch/Mode3_ForBinomialBatch.json")


Model = SC.Circuits.Transmon_Resonators.load(models[ARGS[1]])


ψ  = eval(Meta.parse(ARGS[5]))
ψ = ψ/qt.norm(ψ)

Base.redirect_stdio(stdout = log_file, stderr = log_file) do 
    @info "The args: $ARGS"
    @debug "Got psi"

    c_op_names = []
    if length(ARGS) > 5
        c_op_names = [ARGS[6:end]...]
    end

    c_ops = []
    for name in c_op_names
        @info "Adding $name to C_ops list"
        push!(c_ops, Model.CandD_Ops[name])
    end

    @debug "Got List of C_ops"

    ρ = ψ*ψ'
    start_time = now()

    @info "C and D Ops: $c_op_names"
    @info "Op Sequence: $(Model.Stuff["Drive_Sequences"]["Binomial_Encoding"])"
    other_ds_properties = Dict{Any, Any}("CandD_Ops" => string(c_op_names))


    #Base.redirect_stdio(stdout = log_file2, stderr = log_file2) do 
    SC.Dynamics.RunPulseSequence(Model, ρ, Model.Stuff["Drive_Sequences"]["Binomial_Encoding"], c_ops = c_ops, run_name = run_name, save_path = the_path*"data/", spns = 0.5, other_ds_properties=other_ds_properties)
    #end

    end_time = now()

    @info "Total Run Time: "*string(Dates.canonicalize(end_time - start_time))

    file_name = the_path*"/logs/counter.log" 

    num_done = 0
    if isfile(file_name)
        num_done = open(file_name, "r") do file
            num_done = Meta.parse(read(file, String))
            @info "Opened File, value is $num_done"
            return num_done
        end
        @info "Writing New File with $(num_done+1)"
        open(file_name, "w") do file
            write(file, string(num_done+1))
            close(file)
        end
    else
        open(file_name, "w") do file
            @info "First Run done"
            write(file, string(1))
            close(file)
            @info "Wrote File"
        end
    end
end

