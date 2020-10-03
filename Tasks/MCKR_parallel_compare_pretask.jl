# The MCKR parallel compare pretask.
# This file compares the speed of computation of a MCKR integral computation for the two case of with and without parallelization.

function Do_subtask1(Input_lines)
    # Reading the inputs from the Input_lines and writes a temperory input julia file for the current tasks.
    Temp_input_file=open(joinpath(pwd(),"Temp","MCKR_parallel_compare_input.jl"),"w")
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
    for i=1:length(Input_lines) # Detecting the MC method.
        line=Input_lines[i]
        if length(line)>11 && line[1:11]=="MC_method: "
            write(Temp_input_file,"MC_method=\"$(line[12:end])\"\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    for i=1:length(Input_lines) # Detecting the Sample size.
        line=Input_lines[i]
        if length(line)>13 && line[1:13]=="Sample_size: "
            write(Temp_input_file,"Sample_size=$(line[14:end])\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    for i=1:length(Input_lines) # Detecting number of requested workers.
        line=Input_lines[i]
        if length(line)>16 && line[1:16]=="Workers_number: "
            write(Temp_input_file,"Workers_number=$(line[17:end])\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    close(Temp_input_file)
end

function Do_subtask2()
    # Making a temperory julia file for non-parallel MC computation.
    include(joinpath(pwd(),"Temp","MCKR_parallel_compare_input.jl"))
    Temp_nonparallel_file=open(joinpath(pwd(),"Temp","MCKR_parallel_compare_nonparallel.jl"),"w")
    write(Temp_nonparallel_file,"include(joinpath(pwd(),\"Temp\",\"MCKR_parallel_compare_input.jl\"))
include(joinpath(pwd(),\"Example_MCKR_bank\",MCKR_file_name))
function nonparallel_t(NN)
    st=time_ns()")
    if MC_method=="Simple"
        write(Temp_nonparallel_file,"\n    ans_nonparallel=sumo_with_S(Box,NN)")
    elseif MC_method=="Antithetic"
        write(Temp_nonparallel_file,"\n    ans_nonparallel=sumo_antithetic_with_S(Box,NN)")
    end
    write(Temp_nonparallel_file,"\n    standard_error=ans_nonparallel[2]/NN
    standard_error=standard_error/(NN-1)
    standard_error=sqrt(standard_error)
    st=time_ns()-st
    return(ans_nonparallel[1],standard_error,st/10^9)
end")
    close(Temp_nonparallel_file)
end

function Do_subtask3()
    # Making a temperory julia file for parallel case.
    include(joinpath(pwd(),"Temp","MCKR_parallel_compare_input.jl"))
    Temp_parallel_file=open(joinpath(pwd(),"Temp","MCKR_parallel_compare_parallel.jl"),"w")
    write(Temp_parallel_file,"include(joinpath(pwd(),\"Temp\",\"MCKR_parallel_compare_input.jl\"))
include(joinpath(pwd(),\"Example_MCKR_bank\",MCKR_file_name))
function parallel_t(NN)
    st=time_ns()\n")
    if MC_method=="Simple"
        for i=1:Workers_number
            write(Temp_parallel_file,"    w$(i)ans1=remotecall(sumo_with_S,$(1+i),Box,floor(NN/Workers_number))\n")
        end
    elseif MC_method=="Antithetic"
        for i=1:Workers_number
            write(Temp_parallel_file,"    w$(i)ans1=remotecall(sumo_antithetic_with_S,$(1+i),Box,floor(NN/Workers_number))\n")
        end
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
        st=time_ns()-st
        return(ans2_1,standard_error,st/10^9)
    end")
    close(Temp_parallel_file)
end

function Do_pretask(Input_lines)
    Do_subtask1(Input_lines)
    Do_subtask2()
    Do_subtask3()
end
