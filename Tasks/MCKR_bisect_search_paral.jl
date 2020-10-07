# The MCKR bisect search paral task.
# The MCKR task is going to use the bisect search to find a sub-box entirely inside the region of multistationarity or reject existing of such a box.

using Distributed # A package for parallelization.

# Loading the temp files written by the pretask file.
include(joinpath(pwd(),"Temp","MCKR_bisect_search_paral_input.jl"))
include(joinpath(pwd(),"Example_MCKR_bank",MCKR_file_name))
include(joinpath(pwd(),"Temp","MCKR_bisect_search_paral.jl"))
include(joinpath(pwd(),"Temp","MCKR_integrate_paral.jl"))

function Do_task()
    # The long string variable to recorde the report, step by step.
    Task_output=["\u066d Main computations.\n"]
    push!(Task_output,"\nConsider the bisect search using the $MC_method Monte Carlo integration with N=$Sample_size for the $(MCKR_file_name[1:end-3]) example with the following initial box.\n$Box\nand with the stop criteria as follows;\nsub-box size limit = $Stop_criterion_subbox_size\nlimit on the number of bisecting steps = $Stop_criterion_bisect_step\n")
    # Adding extra workers.
    addprocs(Workers_number)
    # Loading the temporary input file for the extra workers as well.
    @everywhere include(joinpath(pwd(),"Temp","MCKR_bisect_search_paral_input.jl"))
    # Loading the MCKR file for all workers.
    @everywhere include(joinpath(pwd(),"Example_MCKR_bank",MCKR_file_name))
    out1=MC_search(Box,10*Workers_number,1) # Just to warm up Julia about this function.
    out1=MC_search(Box,Sample_size,Stop_criterion_bisect_step)
    if Searching_time_request==1
        push!(Task_output,"\n$(out1[1]).\nIt took $(out1[3]) seconds and computed $(out1[2]) integrals.\nNumber of bisecting step = $(out1[4])\nThe last sub-box = $(out1[5])\nwith I-hat = $(out1[6])\n")
    else
        push!(Task_output,"\n$(out1[1]).\nIt computed $(out1[2]) integrals.\nNumber of bisecting step = $(out1[3])\nThe last sub-box = $(out1[4])\nwith I-hat = $(out1[5])\n")
    end
    if Last_box_info_request==1
        out2=MC_integration(out1[5],10*Workers_number) # just to warm-up Julia for this function.
        out2=MC_integration(out1[5],Sample_size)
        push!(Task_output,"\nComputing the MCKR on the last sub-box again.\nI-hat = $(out2[1]), e-hat = $(out2[2]))\n")
    end
    push!(Task_output,"\n(Parallelization using $Workers_number workers is used)\n\n\u066d End of computations.")
    for i=2:1+Workers_number
        rmprocs(i,waitfor=0)
    end
    return(Task_output)
end
