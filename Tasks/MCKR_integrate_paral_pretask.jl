# The MCKR integrate pretask.
# The MCKR task is going to compute a single Monte Carlo integration for the given example, parameter box and sample size.

function Do_subtask1(Input_lines)
    # Reading the inputs from the Input_lines and writes a temporary input julia file for the current tasks.
    Temp_input_file=open(joinpath(pwd(),"Temp","MCKR_integrate_paral_input.jl"),"w")
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
    # Making a temporary julia file for parallel case.
    include(joinpath(pwd(),"Temp","MCKR_integrate_paral_input.jl"))
    Temp_parallel_file=open(joinpath(pwd(),"Temp","MCKR_integrate_paral.jl"),"w")

    Method_dict=Dict(["Simple",0]=>"sumo",["Simple",1]=>"sumo_with_S",["Antithetic",0]=>"sumo_antithetic",["Antithetic",1]=>"sumo_antithetic_with_S")
    Time_dict=Dict(1=>["    st=time_ns()\n","    st=time_ns()-st\n",",st/10^9"],0=>["","",",0"])
    e_hat_dict=Dict(1=>["[1]","    standard_error=ans2_2/(Workers_number*floor(NN/Workers_number))\n    standard_error=standard_error/((Workers_number*floor(NN/Workers_number))-1)\n    standard_error=sqrt(standard_error)\n",",standard_error"],0=>["","",",0"])
    write(Temp_parallel_file,"include(joinpath(pwd(),\"Temp\",\"MCKR_integrate_paral_input.jl\"))
include(joinpath(pwd(),\"Example_MCKR_bank\",MCKR_file_name))
function MC_integration(B,NN)\n$(Time_dict[integration_time_request][1])")
    for i=1:Workers_number
        write(Temp_parallel_file,"    w$(i)ans1=remotecall($(Method_dict[[MC_method,e_hat_request]]),$(1+i),B,floor(NN/Workers_number))\n")
    end
    for i=1:Workers_number
        write(Temp_parallel_file,"    ans2w$(i)=fetch(w$(i)ans1)\n")
    end
    write(Temp_parallel_file,"    ans2_1=(")
    for i=1:Workers_number-1
        write(Temp_parallel_file,"ans2w$(i)$(e_hat_dict[e_hat_request][1])+")
    end
    write(Temp_parallel_file,"ans2w$Workers_number$(e_hat_dict[e_hat_request][1]))/Workers_number\n")
    if e_hat_request==1
        write(Temp_parallel_file,"    ans2_2=")
        for i=1:Workers_number-1
            write(Temp_parallel_file,"ans2w$(i)[2]+")
        end
        write(Temp_parallel_file,"ans2w$Workers_number[2]\n")
    end
    write(Temp_parallel_file,"$(e_hat_dict[e_hat_request][2])$(Time_dict[integration_time_request][2])    return(ans2_1$(e_hat_dict[e_hat_request][3])$(Time_dict[integration_time_request][3]))
end")
    close(Temp_parallel_file)
end

function Do_pretask(Input_lines)
    Do_subtask1(Input_lines)
    Do_subtask2()
end
