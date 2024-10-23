import QuantumToolbox as qt
using Logging
import CairoMakie as cm
using MiniLoggers
using ProgressMeter
using LoggingExtras
using Revise
using Dates
using YAXArrays


import SuperconductingCavities as SC

df = DateFormat("e-u-d-yy:HH:MM")
t = now()

the_time = string(Dates.format(t, df))


if !isdir("logs/")
    mkdir("logs")
end
log_file = open("logs/CheckingCollapseAndDephasing_with2_gf_sqrt"*the_time*".log", "a")

InfoLogger = MiniLogger(io = log_file, minlevel = MiniLoggers.Info)
#ProgressLogger = MiniLogger(minlevel = LogLevel(-1))
#DebugLogger = MiniLogger(minlevel = MiniLoggers.Debug)

global_logger(InfoLogger)


Base.redirect_stdio(stdout = log_file, stderr = log_file) do 
    Mode3 = SC.Circuits.Transmon_Resonators.load("ModelSaves/Mode3/Mode3.json")

    Ĥ = Mode3.Ĥ
    Ô_D = Mode3.n̂ₜ

    op_params = Dict{Any, Any}("freq_d" => 0, "epsilon" => 0, "Envelope" => "Square", "Envelope Args" => Dict{Any, Any}(), "shift" => 0, "pulse_time" => 2e3)

    ψ = (Mode3.dressed_states[(0,0)]+Mode3.dressed_states[(2,0)])/sqrt(2)

    ρ = ψ*ψ'
    @info Mode3.params
    @info "Running Ramsey"

    SC.Dynamics.RunSingleOperator(Ĥ, Ô_D, ρ, op_params; c_ops = collect(values(Mode3.CandD_Ops)), save_step = true, step_name = "CheckingCollapseAndDephasing_with_2_gf_sqrt", run_name = "CheckingCollapseAndDephasing_with_2_gf_sqrt", op_name = "Ramsey", to_return = "Nothing") 

end