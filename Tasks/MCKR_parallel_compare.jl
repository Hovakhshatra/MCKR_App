# The MCKR parallel compare task.
# This file compares the speed of computation of a MCKR integral computation for the two case of with and without parallelization.

using Distributed # A package for parallelization.
# Loading temperory julia files written by pretask file.
include(joinpath(pwd(),"Temp","MCKR_parallel_compare_input.jl"))
include(joinpath(pwd(),"Temp","MCKR_parallel_compare_nonparallel.jl"))
include(joinpath(pwd(),"Temp","MCKR_parallel_compare_parallel.jl"))

function Do_task()
    # The long string variable to recorde the report, step by step.
    Task_output=["\u066d Information of the computing system.\n"]
    # Adding extra workers.
    addprocs(Workers_number)
    # Loading the temperory input file for the extra workers as well.
    @everywhere include(joinpath(pwd(),"Temp","MCKR_parallel_compare_input.jl"))
    # Loading the MCKR file for all workers.
    @everywhere include(joinpath(pwd(),"Example_MCKR_bank",MCKR_file_name))
    push!(Task_output,"\n\u066d Main computations.\n")
    push!(Task_output,"\nConsider the $MC_method Monte Carlo integration for the $(MCKR_file_name[1:end-3]) example on the box\n$Box\nwith uniform distributions and N=$Sample_size.\n")
    push!(Task_output,"\nWithout parallelization.\n")
    out1_1=nonparallel_t(10) # to warm-up Julia about this function.
    out1_1=nonparallel_t(Sample_size)
    push!(Task_output,"I-hat = $(out1_1[1]), e-hat = $(out1_1[2]), time = $(out1_1[3])\n")
    push!(Task_output,"\nWith parallelization and using $Workers_number workers (workers number 2 till $(1+Workers_number)).\n")
    out2_1=parallel_t(10*Workers_number) # to warm-up Julia about this function.
    out2_1=parallel_t(Sample_size)
    push!(Task_output,"I-hat = $(out2_1[1]), e-hat = $(out2_1[2]), time = $(out2_1[3])\n")
    push!(Task_output,"\n\u066d End of computations.")
    for i=2:1+Workers_number
        rmprocs(i,waitfor=0)
    end
    return(Task_output)
end
