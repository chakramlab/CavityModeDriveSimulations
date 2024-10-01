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

Mode3 = SC.Transmon_Resonators_Loader("ModelSaves/Mode3/Mode3.json");


df = DateFormat("e-u-d-yy:H:M")
t = now()

println(Dates.format(t, df))


c_d_op_combos = [[], ["Mode 3 Collapse"], ["Transmon Collapse"], ["Transmon Dephasing"], ["Mode 3 Collapse", "Transmon Collapse", "Transmon Dephasing"]]

num_combos = length(c_d_op_combos)
for i in 1:length(c_d_op_combos)
    @info "Submitting Job $i/$num_combos"
    c_d_ops_to_use = c_d_op_combos[i]
    
    run_name = "g0_m_e0_Encoding_"*string(Dates.format(t, df));
    for op in c_d_ops_to_use
        run_name = run_name*"_\""*op*"\""
    end

    ψ_string = raw"\"Model.dressed_states[(0,0)] - Model.dressed_states[(1,0)]\""

    string_to_exec = "exec -a $run_name nohup julia BinomialCodeRuns.jl Mode3 $run_name  $ψ_string > Nohup_Outs/$run_name.out &"

    for op in c_d_ops_to_use
        string_to_exec = string_to_exec*" \""*op*"\""
    end

    @info string_to_exec

    run(Cmd(`bash -c $string_to_exec`, ignorestatus=true, detach=true), wait = false)

    if i < num_combos
        min_to_sleep = 5
        time_to_submit_next = Dates.format(now()+Minute(min_to_sleep), "H:M")
        @info "Sleeping for $min_to_sleep, submitting next job at $time_to_submit_next"
        sleep(min_to_sleep*60)
    end
end

@info "Done!"