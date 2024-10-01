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

function tostr(obj)
    io = IOBuffer()
    show(io, "text/plain", obj)
    String(take!(io))
end


Mode3 = Mode3 = SC.Transmon_Resonators_Loader("ModelSaves/Mode3/Mode3.json")

solver_kwargs = Dict{Any, Any}()#"reltol" => 1e-8, "abstol" => 1e-8, "tol"=>1e-8)
ψ = Mode3.dressed_states[(2,0)]
ρ = ψ*ψ'
start_time = now()

with_loss = true
c_ops = []
if with_loss
    c_ops = collect(values(Mode3.CandD_Ops))
end

@info "With Loss: $with_loss"

SC.RunPulseSequence(Mode3, ρ, fill("sb_f0g1", 400), c_ops = c_ops, run_name = "HellaSidebands"*string(with_loss)*"_"*string(now()), solver_kwargs = solver_kwargs, spns = 0.1)
end_time = now()

@info "Total Run Time: "*string(Dates.canonicalize(end_time - start_time))