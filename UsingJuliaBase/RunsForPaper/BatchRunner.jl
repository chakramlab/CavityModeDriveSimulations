import QuantumToolbox as qt
using Logging
import CairoMakie as cm
using MiniLoggers
using Revise
using Dates


import SuperconductingCavities as SC

# MiniLogger(minlevel = MiniLoggers.Info) |> global_logger
InfoLogger = MiniLogger(minlevel = MiniLoggers.Info)
ProgressLogger = MiniLogger(minlevel = LogLevel(-1))
DebugLogger = MiniLogger(minlevel = MiniLoggers.Debug)

global_logger(ProgressLogger)

Mode3 = SC.Circuits.Transmon_Resonators.load("ModelSaves/Mode3/Mode3.json");


df = DateFormat("e-u-d-yy:H:M")
t = now()

println(Dates.format(t, df))


c_d_op_combos = [[], ["Dressed Mode 3 Collapse"], ["Dressed Transmon Collapse"], ["Dressed Transmon Dephasing"], ["Dressed Mode 3 Collapse", "Dressed Transmon Collapse", "Dressed Transmon Dephasing"]]

num_combos = length(c_d_op_combos)

states  = ["g0", "e0", "g0_p_e0", "g0_m_e0", "g0_p_ie0", "g0_m_ie0"]

states_dict = Dict{Any, Any}()
states_dict["g0"] = raw"\"Model.dressed_states[(0,0)]\""
states_dict["g0"] = raw"\"Model.dressed_states[(1,0)]\""
states_dict["g0_p_e0"] = raw"\"Model.dressed_states[(0,0)] + Model.dressed_states[(1,0)]\""
states_dict["g0_m_e0"] = raw"\"Model.dressed_states[(0,0)] - Model.dressed_states[(1,0)]\""
states_dict["g0_p_ie0"] = raw"\"Model.dressed_states[(0,0)] + 1im*Model.dressed_states[(1,0)]\""
states_dict["g0_m_ie0"] = raw"\"Model.dressed_states[(0,0)] - 1im*Model.dressed_states[(1,0)]\""

the_time = string(Dates.format(t, df))
for state in states
    for i in 1:length(c_d_op_combos)
        @info "Submitting Job $i/$num_combos"
        c_d_ops_to_use = c_d_op_combos[i]
        
        run_name = state*"_Encoding_"*the_time;
        for op in c_d_ops_to_use
            run_name = run_name*"_\""*op*"\""
        end

        ψ_string = states_dict[state]

        string_to_exec = "exec -a $run_name nohup julia BinomialCodeRuns.jl Mode3 $run_name  $ψ_string"

        for op in c_d_ops_to_use
            string_to_exec = string_to_exec*" \""*op*"\""
        end

        string_to_exec = string_to_exec*" > Nohup_Outs/$run_name.out &"
        @info string_to_exec

        run(Cmd(`bash -c $string_to_exec`, ignorestatus=true, detach=true), wait = false)

        if i < num_combos
            min_to_sleep = 5
            time_to_submit_next = Dates.format(now()+Minute(min_to_sleep), "H:M")
            @info "Sleeping for $min_to_sleep minutes, submitting next job at $time_to_submit_next"
            sleep(min_to_sleep*60)
        end
    end
    @info "Done submitting jobs for state $state"
    hours_to_sleep = 9
    time_to_submit_next = Dates.format(now()+Hours(hours_to_sleep), "U d H:M")
    @info "Sleeping for $hours_to_sleep hours, moving onto next state at $time_to_submit_next"
    sleep(hours_to_sleep*60*60)
end
@info "Done!"