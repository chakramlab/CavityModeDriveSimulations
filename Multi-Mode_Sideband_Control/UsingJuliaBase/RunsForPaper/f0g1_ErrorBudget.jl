using Pkg
Pkg.activate("SCC", shared=true)

import QuantumToolbox as qt
import SuperconductingCavities as SC
import CairoMakie as cm
using Revise
using Dates
import JSON3

using YAXArrays

using ProgressMeter

import Optim as opt

using Logging
using MiniLoggers

using IJulia
if isdefined(Main, :IJulia)
   Main.IJulia.stdio_bytes[] = 0;
end


# MiniLogger(minlevel = MiniLoggers.Info) |> global_logger

log_file = "logs/f0g1_error_budget.log"

log_file = open(log_file, "a")
InfoLogger = MiniLogger(io = log_file, minlevel = MiniLoggers.Info)
#ProgressLogger = MiniLogger(minlevel = LogLevel(-1))
#DebugLogger = MiniLogger(minlevel = MiniLoggers.Debug)

global_logger(InfoLogger)


df = DateFormat("e-u-d-yy_H_M")
t = now()
name_time = string(Dates.format(t, df))


modes = [6,7,8,9]
Models = [SC.Circuits.Transmon_Resonators.load("ModelSaves/Mode$i/Mode$i.json") for i in modes];

cops_to_do = ["None", "All", "TC", "TD", "CC"]
drive_errors = ["LL", "L", "HL", "None", "HU", "U", "UU"]

Base.redirect_stdio(stdout = log_file, stderr = log_file) do 
    start_time = now()
    for j in 1:length(modes)
        i = modes[j]
        Model = Models[j]
        @info "Doing Mode $(i)"

        f0g1_original = deepcopy(Model.Stuff["op_drive_params"]["sb_f0g1"])
        
        for drive_error in drive_errors
            for cop_to_do in cops_to_do
                @info "Doing loss: $(cop_to_do), drive error: $(drive_error)"
                Model.Stuff["op_drive_params"]["sb_f0g1"] = deepcopy(f0g1_original)
                if drive_error == "LL"
                    Model.Stuff["op_drive_params"]["sb_f0g1"]["epsilon"] *=(1-1/125)
                elseif drive_error == "L"
                    Model.Stuff["op_drive_params"]["sb_f0g1"]["epsilon"] *=(1-1/250)
                elseif drive_error == "HL"
                    Model.Stuff["op_drive_params"]["sb_f0g1"]["epsilon"] *=(1-1/500)
                elseif drive_error == "HU"
                    Model.Stuff["op_drive_params"]["sb_f0g1"]["epsilon"] *=(1+1/500)
                elseif drive_error == "U"
                    Model.Stuff["op_drive_params"]["sb_f0g1"]["epsilon"] *=(1+1/250)
                elseif drive_error == "UU"
                    Model.Stuff["op_drive_params"]["sb_f0g1"]["epsilon"] *=(1+1/125)
                end

                c_ops = []
                if cop_to_do == "TC"
                    c_ops = [Model.CandD_Ops["Transmon Collapse"]]
                elseif cop_to_do == "TD"
                    c_ops = [Model.CandD_Ops["Transmon Dephasing"]]
                elseif cop_to_do == "CC"
                    c_ops = [Model.CandD_Ops["Mode$(i) Collapse"]]
                elseif cop_to_do == "All"
                    c_ops = collect(values(Model.CandD_Ops))
                end

    
                ψ = Model.dressed_states[(2,0)]
                ρ = ψ*ψ'

                other_ds_properties = Dict{Any, Any}()

                #Base.redirect_stdio(stdout = log_file2, stderr = log_file2) do 
                SC.Dynamics.RunPulseSequence(Model, ρ, ["sb_f0g1"], c_ops = c_ops, run_name = "Mode$(i)_Loss_$(cop_to_do)_DE_$(drive_error)", save_path = "Data/f0g1_Error_Budget_$(name_time)/", spns = 2, other_ds_properties=other_ds_properties)
                #end
            end
        end
    end
    end_time = now()
        @info "Total Run Time: "*string(Dates.canonicalize(end_time - start_time))
end
