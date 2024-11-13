import QuantumToolbox as qt
using Logging
using LoggingExtras
import CairoMakie as cm

using Revise
using Dates


import SuperconductingCavities as SC

# MiniLogger(minlevel = MiniLoggers.Info) |> global_logger

df = DateFormat("e-u-d-yy.H.M")
t = now()

the_time = string(Dates.format(t, df))

batch_name = "Binomial_Encoding_"*the_time
the_path = "Data/"*batch_name*"/"

if !isdir(the_path)
    mkdir(the_path)
    mkdir(the_path*"/logs")
    mkdir(the_path*"/data")
end


log_file = the_path*"/logs/BatchRunner.log"
logger = FormatLogger(open(log_file, "a")) do io, args
    df = DateFormat("e-u-d-yy.H.M")
    t = now()
    the_time = string(Dates.format(t, df))
    println(io, the_time, ":", "[", args.level, "] ", args.message )
end

#logger = FileLogger(log_file; append = true)
logger = MinLevelLogger(logger, LogLevel(-1))
global_logger(logger)

Mode3 = SC.Circuits.Transmon_Resonators.load("ModelSaves/Mode3/Mode3.json");

Mode3 = copy(Mode3; name_addon = "ForBinomialBatch")
SC.Utils.save_model(Mode3)

@info "Time of Run: "*String(Dates.format(t, df))


c_d_op_combos = [[], ["Bare Mode 3 Collapse"], ["Bare Transmon Collapse"], ["Bare Transmon Dephasing"], ["Bare Mode 3 Collapse", "Bare Transmon Collapse", "Bare Transmon Dephasing"]]

num_combos = length(c_d_op_combos)

states  = ["g0", "e0", "g0_p_e0", "g0_m_e0", "g0_p_ie0", "g0_m_ie0"]

states_dict = Dict{Any, Any}()
states_dict["g0"] = raw"\"Model.dressed_states[(0,0)]\""
states_dict["e0"] = raw"\"Model.dressed_states[(1,0)]\""
states_dict["g0_p_e0"] = raw"\"Model.dressed_states[(0,0)] + Model.dressed_states[(1,0)]\""
states_dict["g0_m_e0"] = raw"\"Model.dressed_states[(0,0)] - Model.dressed_states[(1,0)]\""
states_dict["g0_p_ie0"] = raw"\"Model.dressed_states[(0,0)] + 1im*Model.dressed_states[(1,0)]\""
states_dict["g0_m_ie0"] = raw"\"Model.dressed_states[(0,0)] - 1im*Model.dressed_states[(1,0)]\""


for state in states
    @info "============================================================================================================================"
    @info "Starting State: $state"
    for i in 1:length(c_d_op_combos)
        @info "Submitting Job $i/$num_combos"
        c_d_ops_to_use = c_d_op_combos[i]
        
        run_name = state#*"_Encoding_"*the_time;
        for op in c_d_ops_to_use
            run_name = run_name*"_\""*op*"\""
        end

        ψ_string = states_dict[state]

        string_to_exec = "exec -a $run_name julia BinomialCodeRuns.jl Mode3 $batch_name $the_path $run_name  $ψ_string"

        ops_used = ""
        for op in c_d_ops_to_use
            ops_used = ops_used*"_"*op
            string_to_exec = string_to_exec*" \""*op*"\""
        end

        string_to_exec = string_to_exec*" &"#*" > $the_path/logs/$(state)$(ops_used).out &"
        @info string_to_exec

        run(Cmd(`bash -c $string_to_exec`, ignorestatus=true, detach=true), wait = false)

        if i < num_combos
            min_to_sleep = 5
            time_to_submit_next = Dates.format(now()+Minute(min_to_sleep), "H:M")
            @info "Sleeping for $min_to_sleep minutes, submitting next job at $time_to_submit_next"
            sleep(5*60)
        end
    end
    @info "Done submitting jobs for state $state. Waiting for them to finish"
    
    file_name = the_path*"/logs/counter.log"
    while true
        num_done = 0
        if isfile(file_name)
            num_done = open(file_name, "r") do file
                num_done = Meta.parse(read(file, String))
                return num_done
            end
        end
        if (num_done % length(c_d_op_combos) == 0) & num_done >0
            @info "All combos done, moving onto next state"
            rm(file_name)
            break
        else
            @info "$num_done completed. Checking again in 10 minutes" 
            sleep(10*60)
        end
    end

end
@info "Done!"