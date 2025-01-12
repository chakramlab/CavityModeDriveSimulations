using Pkg
Pkg.activate("SCC", shared=true)

import QuantumToolbox as qt
using Logging
using LoggingExtras
using DimensionalData
using YAXArrays
using Dates

using MiniLoggers

import SuperconductingCavities as SC

# MiniLogger(minlevel = MiniLoggers.Info) |> global_logger

log_file = "logs/FloquetSweepingForIonization.log"

log_file = open(log_file, "a")
InfoLogger = MiniLogger(io = log_file, minlevel = MiniLoggers.Info)
#ProgressLogger = MiniLogger(minlevel = LogLevel(-1))
#DebugLogger = MiniLogger(minlevel = MiniLoggers.Debug)

global_logger(InfoLogger)


Model = SC.Circuits.Transmon_Resonators.load("ModelSaves/Mode3/Mode3.json")

df = DateFormat("e-u-d-yy_HH_MM")
t = now()
name_time = string(Dates.format(t, df))

Base.redirect_stdio(stdout = log_file, stderr = log_file) do 
    freq_d = Model.dressed_energies[(0,1)]-Model.dressed_energies[(2,0)]
    #εs = collect(LinRange(0.7, 1.1, 200));
    #shifts = collect(LinRange(0.03, 0.1, 200))
    εs = collect(LinRange(0.7, 1.1, 2));
    shifts = collect(LinRange(0.03, 0.1, 11))
    states = keys(Model.dressed_energies)
    quasienergies = fill(0.0, Dim{:state}(string.(states)), Dim{:eps}(εs), Dim{:shift}(shifts))

        
    start_time = now()

    for j in 1:length(shifts)
        step_start = now()
        @info "Doing Shift Number: $(j)"
        params = []
        for i in 1:length(εs)
            push!(params, Dict{Any,Any}())
            params[end]["ε"] = εs[i]
            params[end]["ν"] = shifts[j]+freq_d
        end
    
        floq_sweep_res = SC.Dynamics.Floquet_t0_Sweep(Model, params; states_to_track=Model.dressed_states);
    
        
    
        for state in states
            #quasienergies[shifts[j]][state] = [floq_sweep_res[State = At(string(state))].data[i]["Quasienergies"] for i in 1:length(εs)]
            floq_dat = [floq_sweep_res[State = At(string(state))].data[i]["Quasienergies"] for i in 1:length(εs)]
            @info floq_dat
            quasienergies[state = At(string(state)), shift=At(shifts[j])] = floq_dat
        end
        step_end = now()
        @info "Step Run Time: "*string(Dates.canonicalize(step_end - step_start))
        
        avg_time = Dates.Millisecond(round(Int, Dates.value(step_end - start_time)/j))
        expected_end = avg_time*(length(shifts)-j)
        @info "Average Time Per Step: "*string(Dates.canonicalize(avg_time))
        @info "Finishing in: "*string(Dates.canonicalize(expected_end))
        @info "--------------------------------------------------------------------------------"
        
    
    end

    end_time = now()

    @info "Total Run Time: "*string(Dates.canonicalize(end_time - start_time))

    savecube(YAXArray(quasienergies.dims, quasienergies.data), "Data/Quasienergies_"*name_time*".nc", overwrite=true)
    
end
