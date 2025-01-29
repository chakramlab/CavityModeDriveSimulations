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

df = DateFormat("e-u-d-yy-HH-MM")
t = now()

the_time = string(Dates.format(t, df))

log_file = open("logs/CalibratingRamps_"*the_time*".log", "a")

InfoLogger = MiniLogger(io = log_file, minlevel = MiniLoggers.Info)
#ProgressLogger = MiniLogger(minlevel = LogLevel(-1))
#DebugLogger = MiniLogger(minlevel = MiniLoggers.Debug)

global_logger(InfoLogger)



Base.redirect_stdio(stdout = log_file, stderr = log_file) do 

    Mode3 = SC.Circuits.Transmon_Resonators.load("ModelSaves/Mode3_ManyRamps/Mode3_ManyRamps.json")

    Ramp_Times = sort(vcat(1.0./collect(1:0.5:10), collect(1:1:20), collect(20:10:100)))
    envelopes = ["Bump_Ramp", "Sine_Squared_Ramp", "Gaussian_Ramp"]

    state1 = (2,0)
    state2 = (0,1)

    ψ1 = Mode3.dressed_states[state1]
    ψ2 = Mode3.dressed_states[state2]

    freq_d = Mode3.dressed_energies[state2]-Mode3.dressed_energies[state1]

    ε = 0.78
    stark_shift = 0.042200385158
    
    SC.Utils.save_model(copy(Mode3; name_addon = "Backup_"*the_time))

    for envelope in envelopes
        for rt in Ramp_Times
            t_range = [180, 180+2*rt]
            for chirp in [false]#, false]
                @info "Doing $envelope with rt: $rt, chirp is $chirp"
                drive_name = "f0g1"
                if envelope == "Sine_Squared_Ramp"
                    drive_name = drive_name*"_SS_"*string(rt)
                end
                if envelope == "Bump_Ramp"
                    drive_name = drive_name*"_B_"*string(rt)
                end
                if envelope == "Gaussian_Ramp"
                    drive_name = drive_name*"_G_"*string(rt)
                end

                if chirp
                    drive_name = drive_name*"_chirped"
                end
                
                envelope_args = Dict{Any, Any}("ramp_time" => rt, "pulse_time" => 0)

                chirp_params = nothing
                if chirp
                    chirp_params = Mode3.Stuff["Chirp Params"]["fn_gn+1"][1]["param"]
                end

                Mode3.Stuff["op_drive_params"][drive_name] = SC.Dynamics.OptimizePulse(Mode3, ψ1, ψ2, ε, freq_d, stark_shift, t_range, envelope, envelope_args, levels = 6, samples_per_level = 7, chirp_params = chirp_params)
                @info ""
                @info "===================================================================================================================================================================================="        
                @info "===================================================================================================================================================================================="        
                @info ""
                @info ""
                SC.Utils.save_model(Mode3)
            end
        end
    end
    @info "DONE :D"
end