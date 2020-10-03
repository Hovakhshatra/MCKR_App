# The MCKR method comparirion table parallelized pretask file.
# This task collect the I-hat, e-hat and integration time for the Simple and Antithetic MC integration for the sample sizes starting from Start_sample_size to End_sample_size by a geometrix Step_sample_size.


function Do_subtask1(Input_lines)
    # Reading the inputs from the Input_lines and writes a temperory input julia file for the current tasks.
    Temp_input_file=open(joinpath(pwd(),"Temp","MCKR_method_comparison_table_paral_input.jl"),"w")
    for i=1:length(Input_lines) # Detecting the MCKR file.
        line=Input_lines[i]
        if length(line)>16 && line[1:16]=="MCKR_file_name: "
            write(Temp_input_file,"MCKR_file_name=\"$(line[17:end]).jl\"\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    for i=1:length(Input_lines) # Detecting the parameter box.
        line=Input_lines[i]
        if length(line)>5 && line[1:5]=="Box: "
            write(Temp_input_file,"Box=[$(line[6:end])]\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    for i=1:length(Input_lines) # Detecting the Sample size.
        line=Input_lines[i]
        if length(line)>13 && line[1:13]=="Sample_size: "
            write(Temp_input_file,"Sample_size_range=[$(line[14:end])]\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    for i=1:length(Input_lines)
        line=Input_lines[i]
        if length(line)>16 && line[1:16]=="Workers_number: " # Detecting number of requested workers.
            if line[17:end]=="cpu_number"
                write(Temp_input_file,"Workers_number=$(length(Sys.cpu_info()))\n")
            elseif line[17:end]=="half_cpu_number"
                write(Temp_input_file,"Workers_number=$(Int(floor(length(Sys.cpu_info())/2)))\n")
            else
                write(Temp_input_file,"Workers_number=$(line[17:end])\n")
            end
            deleteat!(Input_lines,i)
            break
        end
    end
    close(Temp_input_file)
end

function Do_subtask2()
    # Making a temperory julia file for single parallel integration with e-hat.
    include(joinpath(pwd(),"Temp","MCKR_method_comparison_table_paral_input.jl"))
    Temp_parallel_file=open(joinpath(pwd(),"Temp","MCKR_integrate_paral.jl"),"w")
    write(Temp_parallel_file,"include(joinpath(pwd(),\"Temp\",\"MCKR_method_comparison_table_paral_input.jl\"))
include(joinpath(pwd(),\"Example_MCKR_bank\",MCKR_file_name))
function MC_simple(B,NN)
    st1=time_ns()\n")
    for i=1:Workers_number
        write(Temp_parallel_file,"    w$(i)ans1=remotecall(sumo_with_S,$(1+i),B,floor(NN/Workers_number))\n")
    end
    for i=1:Workers_number
        write(Temp_parallel_file,"    ans2w$(i)=fetch(w$(i)ans1)\n")
    end
    write(Temp_parallel_file,"    ans2_1=(")
    for i=1:Workers_number-1
        write(Temp_parallel_file,"ans2w$(i)[1]+")
    end
    write(Temp_parallel_file,"ans2w$Workers_number[1])/Workers_number\n")
    write(Temp_parallel_file,"    ans2_2=")
    for i=1:Workers_number-1
        write(Temp_parallel_file,"ans2w$(i)[2]+")
    end
    write(Temp_parallel_file,"ans2w$Workers_number[2]
    standard_error=ans2_2/(Workers_number*floor(NN/Workers_number))
    standard_error=standard_error/((Workers_number*floor(NN/Workers_number))-1)
    standard_error=sqrt(standard_error)
    st2=time_ns()
    return(ans2_1,standard_error,(st2-st1)/10^9)
end\n")
    write(Temp_parallel_file,"function MC_antithetic(B,NN)
    st1=time_ns()\n")
        for i=1:Workers_number
            write(Temp_parallel_file,"    w$(i)ans1=remotecall(sumo_antithetic_with_S,$(1+i),B,floor(NN/Workers_number))\n")
        end
        for i=1:Workers_number
            write(Temp_parallel_file,"    ans2w$(i)=fetch(w$(i)ans1)\n")
        end
        write(Temp_parallel_file,"    ans2_1=(")
        for i=1:Workers_number-1
            write(Temp_parallel_file,"ans2w$(i)[1]+")
        end
        write(Temp_parallel_file,"ans2w$Workers_number[1])/Workers_number\n")
        write(Temp_parallel_file,"    ans2_2=")
        for i=1:Workers_number-1
            write(Temp_parallel_file,"ans2w$(i)[2]+")
        end
        write(Temp_parallel_file,"ans2w$Workers_number[2]
    standard_error=ans2_2/(Workers_number*floor(NN/Workers_number))
    standard_error=standard_error/((Workers_number*floor(NN/Workers_number))-1)
    standard_error=sqrt(standard_error)
    st2=time_ns()
    return(ans2_1,standard_error,(st2-st1)/10^9)
end")
    close(Temp_parallel_file)
end

function Do_pretask(Input_lines)
    Do_subtask1(Input_lines)
    Do_subtask2()
end
