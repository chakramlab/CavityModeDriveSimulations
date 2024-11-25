import QuantumToolbox as qt
using Logging
using LoggingExtras

using Dates

using MiniLoggers

import SuperconductingCavities as SC

# MiniLogger(minlevel = MiniLoggers.Info) |> global_logger

log_file = "logs/f0g1_fidelities_upper.log"

log_file = open(log_file, "a")
InfoLogger = MiniLogger(io = log_file, minlevel = MiniLoggers.Info)
#ProgressLogger = MiniLogger(minlevel = LogLevel(-1))
#DebugLogger = MiniLogger(minlevel = MiniLoggers.Debug)

global_logger(InfoLogger)


Models = [SC.Circuits.Transmon_Resonators.load("ModelSaves/Mode$i/Mode$i.json") for i in 1:10];


Base.redirect_stdio(stdout = log_file, stderr = log_file) do 
    for i in 1:10
        Model = Models[i]

        for key in keys(Model.Stuff["op_drive_params"])
            if string(key[1]) == "s"
                @info "Changing epsilon for $key"
                Model.Stuff["op_drive_params"][key]["epsilon"]*= (1+1/250)
            end
        end

        c_ops = collect(values(Model.CandD_Ops))
        ψ = Model.dressed_states[(2,0)]
        ρ = ψ*ψ'
        start_time = now()

        other_ds_properties = Dict{Any, Any}()


        #Base.redirect_stdio(stdout = log_file2, stderr = log_file2) do 
        SC.Dynamics.RunPulseSequence(Model, ρ, ["sb_f0g1"], c_ops = c_ops, run_name = "f0g1_fidelities_mode_$(i)_upper", save_path = "Data/", spns = 2, other_ds_properties=other_ds_properties)
        #end

        end_time = now()

        @info "Total Run Time: "*string(Dates.canonicalize(end_time - start_time))
    end
end
