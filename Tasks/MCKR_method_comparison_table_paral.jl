# The MCKR method comparirion table parallelized task file.
# This task collect the I-hat, e-hat and integration time for the Simple and Antithetic MC integration for the sample sizes starting from Start_sample_size to End_sample_size by a geometrix Step_sample_size.

using Distributed # A package for parallelization.

# Loading the temp files written by the pretask file.
include(joinpath(pwd(),"Temp","MCKR_method_comparison_table_paral_input.jl"))
include(joinpath(pwd(),"Temp","MCKR_integrate_paral.jl"))

Sample_size_start=Sample_size_range[1]
while Sample_size_start<=Sample_size_range[3]
    if Sample_size_start>2*Workers_number
        break
    end
    global Sample_size_start=Sample_size_start*Sample_size_range[2]
end
Sample_size_step_size=Sample_size_range[2]
Sample_size_step_number=Int(floor(log(Sample_size_step_size,(Sample_size_range[3]+1)/Sample_size_start)))+1

function Do_task()
    # The long string variable to recorde the report, step by step.
    Task_output=["\u066d Main computations.\n"]
    push!(Task_output,"\nHere we compare e-hat and the computation time of the simple and antithetic Monte Carlo integrations for the $(MCKR_file_name[1:end-3]) example on the following box\n$Box\nThe considered sample size runs from $Sample_size_start to $(Sample_size_range[3]) with $(Sample_size_range[2])-fold increases.\n\n")
    # Adding extra workers.
    addprocs(Workers_number)
    # Loading the temporary input file for the extra workers as well.
    @everywhere include(joinpath(pwd(),"Temp","MCKR_method_comparison_table_paral_input.jl"))
    # Loading the MCKR file for all workers.
    @everywhere include(joinpath(pwd(),"Example_MCKR_bank",MCKR_file_name))
    out1=MC_simple(Box,10*Workers_number) # Just to warm up Julia about this function.
    out2=MC_antithetic(Box,10*Workers_number) # Just to warm up Julia about this function.
    for ii=1:Sample_size_step_number
        Sample_size_current=Sample_size_start*Sample_size_step_size^(ii-1)
        push!(Task_output,"N = $Sample_size_current\n")
        print("N = $Sample_size_current\n")
        out1=MC_simple(Box,Sample_size_start*Sample_size_step_size^(ii-1))
        push!(Task_output," With Simple Monte Carlo:\n  $(out1[1]) = I-hat\n  $(out1[2]) = e-hat\n  $(out1[3]) = computation time.\n")
        print(" With Simple Monte Carlo:\n  $(out1[1]) = I-hat\n  $(out1[2]) = e-hat\n  $(out1[3]) = computation time.\n")
        out2=MC_antithetic(Box,Sample_size_start*Sample_size_step_size^(ii-1))
        push!(Task_output," With Antithetic Monte Carlo:\n  $(out2[1]) = I-hat\n  $(out2[2]) = e-hat\n  $(out2[3]) = computation time.\n")
        print(" With Antithetic Monte Carlo:\n  $(out2[1]) = I-hat\n  $(out2[2]) = e-hat\n  $(out2[3]) = computation time.\n")
        if out2[3]*out2[2]!=0
            push!(Task_output,"   The efficiency of the antithetic to the simple method = $((out1[3]*out1[2])/(out2[3]*out2[2]))\n")
            print("   The efficiency of the antithetic to the simple method = $((out1[3]*out1[2])/(out2[3]*out2[2]))\n")
        end
    end
    for i=2:1+Workers_number
        rmprocs(i,waitfor=0)
    end
    push!(Task_output,"\n\u066d End of computations.")
    return(Task_output)
end
