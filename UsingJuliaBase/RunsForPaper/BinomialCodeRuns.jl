import QuantumToolbox as qt
using Logging

using MiniLoggers
using Dates


import SuperconductingCavities as SC

# MiniLogger(minlevel = MiniLoggers.Info) |> global_logger
InfoLogger = MiniLogger(minlevel = MiniLoggers.Info)
ProgressLogger = MiniLogger(minlevel = LogLevel(-1))
DebugLogger = MiniLogger(minlevel = MiniLoggers.Debug)

global_logger(DebugLogger)


println(ARGS)


models = Dict{Any, Any}("Mode3" => "ModelSaves/Mode3/Mode3.json")


Model = SC.Circuits.Transmon_Resonators.load(models[ARGS[1]])

@debug "Loaded Model"

run_name = ARGS[2]
ψ = eval(Meta.parse(ARGS[3]))
ψ = ψ/qt.norm(ψ)

@debug "Got psi"

c_op_names = []
if length(ARGS) > 3
    c_op_names = [ARGS[4:end]...]
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

SC.Dynamics.RunPulseSequence(Model, ρ, Model.Stuff["Drive_Sequences"]["Binomial_Encoding"], c_ops = c_ops, run_name = run_name, spns = 0.5, other_ds_properties=other_ds_properties)
end_time = now()

@info "Total Run Time: "*string(Dates.canonicalize(end_time - start_time))