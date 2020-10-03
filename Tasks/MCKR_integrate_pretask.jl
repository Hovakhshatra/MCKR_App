# The MCKR integrate pretask.
# The MCKR task is going to compute a single Monte Carlo integration for the given example, parameter box and sample size.

function Do_subtask1(Input_lines)
    # Reading the inputs from the Input_lines and writes a temperory input julia file for the current tasks.
    Temp_input_file=open(joinpath(pwd(),"Temp","MCKR_integrate_input.jl"),"w")
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
    e_hat_request_check=0 # Checking if the standard error is requested.
    for i=1:length(Input_lines)
        line=Input_lines[i]
        if length(line)==22 && line[1:22]=="Include standard error"
            e_hat_request_check=1
            write(Temp_input_file,"e_hat_request=1\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    if e_hat_request_check==0
        write(Temp_input_file,"e_hat_request=0\n")
    end
    integration_time_request_check=0 # Checking if the standard error is requested.
    for i=1:length(Input_lines)
        line=Input_lines[i]
        if length(line)==24 && line[1:24]=="Include integration time"
            integration_time_request_check=1
            write(Temp_input_file,"integration_time_request=1\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    if integration_time_request_check==0
        write(Temp_input_file,"integration_time_request=0\n")
    end
    close(Temp_input_file)
end

function Do_subtask2()
    # Making a temperory julia file for non-parallel MC computation.
    include(joinpath(pwd(),"Temp","MCKR_integrate_input.jl"))
    Temp_nonparallel_file=open(joinpath(pwd(),"Temp","MCKR_integrate.jl"),"w")
    write(Temp_nonparallel_file,"include(joinpath(pwd(),\"Temp\",\"MCKR_integrate_input.jl\"))
include(joinpath(pwd(),\"Example_MCKR_bank\",MCKR_file_name))
function MC_integration(B,NN)")
    if integration_time_request==1
    write(Temp_nonparallel_file,"    st=time_ns()")
        if e_hat_request==1
            if MC_method=="Simple"
                write(Temp_nonparallel_file,"\n    ans_nonparallel=sumo_with_S(B,NN)")
            elseif MC_method=="Antithetic"
                write(Temp_nonparallel_file,"\n    ans_nonparallel=sumo_antithetic_with_S(B,NN)")
            end
            write(Temp_nonparallel_file,"\n    standard_error=ans_nonparallel[2]/NN
    standard_error=standard_error/(NN-1)
    standard_error=sqrt(standard_error)
    st=time_ns()-st
    return(ans_nonparallel[1],standard_error,st/10^9)
end")
        else
            if MC_method=="Simple"
                write(Temp_nonparallel_file,"\n    ans_nonparallel=sumo(B,NN)")
            elseif MC_method=="Antithetic"
                write(Temp_nonparallel_file,"\n    ans_nonparallel=sumo_antithetic(B,NN)")
            end
            write(Temp_nonparallel_file,"\n    st=time_ns()-st
    return(ans_nonparallel[1],st/10^9)
end")
        end
    else
        if e_hat_request==1
            if MC_method=="Simple"
                write(Temp_nonparallel_file,"\n    ans_nonparallel=sumo_with_S(B,NN)")
            elseif MC_method=="Antithetic"
                write(Temp_nonparallel_file,"\n    ans_nonparallel=sumo_antithetic_with_S(B,NN)")
            end
            write(Temp_nonparallel_file,"\n    standard_error=ans_nonparallel[2]/NN
    standard_error=standard_error/(NN-1)
    standard_error=sqrt(standard_error)
    return(ans_nonparallel[1],standard_error)
end")
        else
            if MC_method=="Simple"
                write(Temp_nonparallel_file,"\n    ans_nonparallel=sumo(B,NN)")
            elseif MC_method=="Antithetic"
                write(Temp_nonparallel_file,"\n    ans_nonparallel=sumo_antithetic(B,NN)")
            end
            write(Temp_nonparallel_file,"\n    return(ans_nonparallel[1])
end")
        end
    end
    close(Temp_nonparallel_file)
end

function Do_pretask(Input_lines)
    Do_subtask1(Input_lines)
    Do_subtask2()
end
