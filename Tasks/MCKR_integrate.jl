# The MCKR integrate task (without parallelization).
# The MCKR task is going to compute a single Monte Carlo integration for the given example, parameter box and sample size.

# The temp files are written by the pretask file.
include(joinpath(pwd(),"Temp","MCKR_integrate_input.jl"))
include(joinpath(pwd(),"Temp","MCKR_integrate.jl"))

function Do_task()
    # The long string variable to recorde the report, step by step.
    Task_output=["\u066d Main computations.\n"]
    push!(Task_output,"\nConsider the $MC_method Monte Carlo integration for the $(MCKR_file_name[1:end-3]) example on the box\n$Box\nand N=$Sample_size.\n")
    out1=MC_integration(Box,100) # Just to warm up Julia about this function.
    out1=MC_integration(Box,Sample_size)
    if integration_time_request==1
        if e_hat_request==1
            push!(Task_output,"\nI-hat = $(out1[1]), e-hat = $(out1[2]), integration time = $(out1[3]) seconds\n")
        else
            push!(Task_output,"\nI-hat = $(out1[1]), integration time = $(out1[3]) seconds\n")
        end
    else
        if e_hat_request==1
            push!(Task_output,"\nI-hat = $(out1[1]), e-hat = $(out1[2])\n")
        else
            push!(Task_output,"\nI-hat = $(out1[1])\n")
        end
    end
    push!(Task_output,"\n\u066d End of computations.")
    return(Task_output)
end
