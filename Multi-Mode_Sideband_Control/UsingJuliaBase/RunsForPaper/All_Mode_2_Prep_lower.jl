import QuantumToolbox as qt
using Logging
using LoggingExtras

using Dates

using MiniLoggers

import SuperconductingCavities as SC

# MiniLogger(minlevel = MiniLoggers.Info) |> global_logger

log_file = "logs/All_Mode_2_Prep_lower.log"

log_file = open(log_file, "a")
InfoLogger = MiniLogger(io = log_file, minlevel = MiniLoggers.Info)
#ProgressLogger = MiniLogger(minlevel = LogLevel(-1))
#DebugLogger = MiniLogger(minlevel = MiniLoggers.Debug)

global_logger(InfoLogger)



df = DateFormat("e-u-d-yy_H_M")
t = now()
name_time = string(Dates.format(t, df))

Base.redirect_stdio(stdout = log_file, stderr = log_file) do 
    for i in 1:10
        Model = SC.Circuits.Transmon_Resonators.load("ModelSaves/Mode$(i)/Mode$(i).json");

        for key in keys(Model.Stuff["op_drive_params"])
            if string(key[1]) == "s"
                @info "Changing epsilon for $key"
                Model.Stuff["op_drive_params"][key]["epsilon"]*= (1-1/250)
            end
        end
        
        c_ops = collect(values(Model.CandD_Ops))
        ψ = Model.dressed_states[(0,0)]
        ρ = ψ*ψ'
        start_time = now()

        other_ds_properties = Dict{Any, Any}()


        #Base.redirect_stdio(stdout = log_file2, stderr = log_file2) do 
        SC.Dynamics.RunPulseSequence(Model, ρ, Model.Stuff["Drive_Sequences"]["Prep_2"], c_ops = c_ops, run_name = "Mode$(i)_Prep_2_lower_$name_time", save_path = "Data/", spns = 2, other_ds_properties=other_ds_properties)
        #end

        end_time = now()

        @info "Total Run Time: "*string(Dates.canonicalize(end_time - start_time))
    end
end