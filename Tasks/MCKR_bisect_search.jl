# The MCKR bisect search pretask.
# The MCKR task is going to use the bisect search to find a sub-box entirely inside the region of multistationarity or reject existing of such a box.

# Loading the temp files written by the pretask file.
include(joinpath(pwd(),"Temp","MCKR_bisect_search_input.jl"))
include(joinpath(pwd(),"Example_MCKR_bank",MCKR_file_name))
include(joinpath(pwd(),"Temp","MCKR_bisect_search.jl"))

function Do_task()
    # The long string variable to recorde the report, step by step.
    Task_output=["\u066d Main computations.\n"]
    push!(Task_output,"\nConsider the bisect search using the $MC_method Monte Carlo integration with N=$Sample_size for the $(MCKR_file_name[1:end-3]) example with the following initial box.\n$Box\nand with the stop criteria as follows;\nsub-box size limit = $Stop_criterion_subbox_size\nlimit on the number of bisecting steps = $Stop_criterion_bisect_step\n")
    out1=MC_search(Box,Sample_size)
    if Searching_time_request==1
        push!(Task_output,"\n$(out1[1]).\nIt took $(out1[3]) seconds and computed $(out1[2]) integrals.\nNumber of bisecting step = $(out1[4])\nThe last sub-box = $(out1[5])\nwith I-hat = $(out1[6])\n")
    else
        push!(Task_output,"\n$(out1[1]).\nIt computed $(out1[2]) integrals.\nNumber of bisecting step = $(out1[3])\nThe last sub-box = $(out1[4])\nwith I-hat = $(out1[5])\n")
    end
    if Last_box_info_request==1
        if MC_method=="Antithetic"
            out2=sumo_antithetic_with_S(out1[5],100) # just to warm-up Julia for this function.
            out2=sumo_antithetic_with_S(out1[5],Sample_size)
        elseif MC_method=="Simple"
            out2=sumo_with_S(out1[5],100) # just to warm-up Julia for this function.
            out2=sumo_with_S(out1[5],Sample_size)
        end
        Standard_error=out2[2]/Sample_size
        Standard_error=Standard_error/(Sample_size-1)
        Standard_error=sqrt(Standard_error)
        push!(Task_output,"\nComputing the MCKR on the last sub-box again.\nI-hat = $(out2[1]), e-hat = $Standard_error\n")
    end
    push!(Task_output,"\n\u066d End of computations.")
    return(Task_output)
end
