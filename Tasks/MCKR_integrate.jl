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
    Output_dict=Dict([1,1]=>"\nI-hat = $(out1[1]), e-hat = $(out1[2]), integration time = $(out1[3]) seconds\n",
                     [1,0]=>"\nI-hat = $(out1[1]), integration time = $(out1[3]) seconds\n",
                     [0,1]=>"\nI-hat = $(out1[1]), e-hat = $(out1[2])\n",
                     [0,0]=>"\nI-hat = $(out1[1])\n")
    push!(Task_output,Output_dict[[integration_time_request,e_hat_request]])
    push!(Task_output,"\n\u066d End of computations.")
    return(Task_output)
end
