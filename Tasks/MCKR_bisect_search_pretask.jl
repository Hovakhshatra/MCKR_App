# The MCKR bisect search pretask.
# The MCKR task is going to use the bisect search to find a sub-box entirely inside the region of multistationarity or reject existing of such a box.

function Do_subtask1(Input_lines)
    # Reading the inputs from the Input_lines and writes a temperory input julia file for the current tasks.
    Temp_input_file=open(joinpath(pwd(),"Temp","MCKR_bisect_search_input.jl"),"w")
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
    for i=1:length(Input_lines) # Detecting the stop criterion with respect to the sob-box size.
        line=Input_lines[i]
        if length(line)>28 && line[1:28]=="Stop_criterion_subbox_size: "
            write(Temp_input_file,"Stop_criterion_subbox_size=[$(line[29:end])]\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    for i=1:length(Input_lines) # Detecting the stop criterion with respect to the number of bisect steps.
        line=Input_lines[i]
        if length(line)>28 && line[1:28]=="Stop_criterion_bisect_step: "
            write(Temp_input_file,"Stop_criterion_bisect_step=$(line[29:end])\n")
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
    Last_box_info_request_check=0 # Checking if the I-hat and e-hat of the last box is asked.
    for i=1:length(Input_lines)
        line=Input_lines[i]
        if length(line)==33 && line[1:33]=="Include last box integration info"
            Last_box_info_request_check=1
            write(Temp_input_file,"Last_box_info_request=1\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    if Last_box_info_request_check==0
        write(Temp_input_file,"Last_box_info_request=0\n")
    end
    Searching_time_request_check=0 # Checking if the computation time of the bisect search is asked.
    for i=1:length(Input_lines)
        line=Input_lines[i]
        if length(line)==22 && line[1:22]=="Include searching time"
            Searching_time_request_check=1
            write(Temp_input_file,"Searching_time_request=1\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    if Searching_time_request_check==0
        write(Temp_input_file,"Searching_time_request=0\n")
    end
    close(Temp_input_file)
end

function Do_subtask2()
    # Making a temperory julia file for bisect search.
    include(joinpath(pwd(),"Temp","MCKR_bisect_search_input.jl"))
    Temp_bisect_search_file=open(joinpath(pwd(),"Temp","MCKR_bisect_search.jl"),"w")
    write(Temp_bisect_search_file,"include(joinpath(pwd(),\"Temp\",\"MCKR_bisect_search_input.jl\"))
include(joinpath(pwd(),\"Example_MCKR_bank\",MCKR_file_name))
function SameBox(B) # To copy a box without dependency between the copy and the original.
    return(deepcopy(B))
end\n")
    write(Temp_bisect_search_file,"function MC_search(B,NN)\n")
    if MC_method=="Simple"
        write(Temp_bisect_search_file,"    sumo(B,100) # Just to warn up Julia about ths function\n")
    elseif MC_method=="Antithetic"
        write(Temp_bisect_search_file,"    sumo_antithetic(B,100) # Just to warn up Julia about ths function\n")
    end
    if Searching_time_request==1
        write(Temp_bisect_search_file,"    st1=time_ns()
    Box_dimension = length(B)
    Index_list=[]
    for i=1:Box_dimension
        push!(Index_list,i)
    end
    for i=1:Box_dimension-1
        push!(Index_list,i)
    end
    Index1=0
    Index2=0\n")
        if MC_method=="Simple"
            write(Temp_bisect_search_file,"    ans0 = sumo(B,NN)\n")
        elseif MC_method=="Antithetic"
            write(Temp_bisect_search_file,"    ans0=sumo_antithetic(B,NN)\n")
        end
        write(Temp_bisect_search_file,"    integral_number=1
    if ans0 < 1.05
        return \"No multistationary sub-box\",integral_number,(time_ns()-st1)/(10^9),Index2,B,ans0
    elseif ans0 > 2.95
        return \"Found multistationary sub-box\",integral_number,(time_ns()-st1)/(10^9),Index2,B,ans0
    end
    @inbounds while Index2 != Stop_criterion_bisect_step
        if Index1 == Box_dimension
            Index1 = 1
        else
            Index1 += 1
        end
        Size_stop=1
        for i=Index1:Index1+Box_dimension-1
            if B[Index_list[i]][2]-B[Index_list[i]][1] > Stop_criterion_subbox_size[Index_list[i]]
                Size_stop=0
                Index2 += 1
                B1 = SameBox(B)
                B1[Index_list[i]][2] = (B[Index_list[i]][1] + B[Index_list[i]][2]) / 2\n")
        if MC_method=="Simple"
            write(Temp_bisect_search_file,"                ans1 = sumo(B1,NN)\n")
        elseif MC_method=="Antithetic"
            write(Temp_bisect_search_file,"                ans1=sumo_antithetic(B1,NN)\n")
        end
        write(Temp_bisect_search_file,"                B2 = SameBox(B)
                B2[Index_list[i]][1] = (B[Index_list[i]][1] + B[Index_list[i]][2]) / 2\n")
        if MC_method=="Simple"
            write(Temp_bisect_search_file,"                ans2 = sumo(B2,NN)\n")
        elseif MC_method=="Antithetic"
            write(Temp_bisect_search_file,"                ans2=sumo_antithetic(B2,NN)\n")
        end
        write(Temp_bisect_search_file,"                integral_number+=2
                if ans1 > ans2
                    B=B1
                    ans0=ans1
                else
                    B=B2
                    ans0=ans2
                end
                if ans0 < 1.05
                    return \"No multistationary sub-box\",integral_number,(time_ns()-st1)/(10^9),Index2,B,ans0
                elseif ans0 > 2.95
                    return \"Found multistationary sub-box\",integral_number,(time_ns()-st1)/(10^9),Index2,B,ans0
                end
                break
            end
        end
        if Size_stop==1
            return \"Termination condition on the size of sub-boxes reached.\",integral_number,(time_ns()-st1)/(10^9),Index2,B,ans0
        end
    end
    st2=time_ns()
    return \"Termination condition on the number of bisect steps reached.\",integral_number,(time_ns()-st1)/(10^9),Index2,B,ans0
end\n")
    else
        write(Temp_bisect_search_file,"    Box_dimension = length(B)
    Index_list=[]
    for i=1:Box_dimension
        push!(Index_list,i)
    end
    for i=1:Box_dimension-1
        push!(Index_list,i)
    end
    Index1=0
    Index2=0\n")
        if MC_method=="Simple"
            write(Temp_bisect_search_file,"    ans0 = sumo(B,NN)\n")
        elseif MC_method=="Antithetic"
            write(Temp_bisect_search_file,"    ans0=sumo_antithetic(B,NN)\n")
        end
        write(Temp_bisect_search_file,"    integral_number=1
    if ans0 < 1.05
        return \"No multistationary sub-box\",integral_number,Index2,B,ans0
    elseif ans0 > 2.95
        return \"Found multistationary sub-box\",integral_number,Index2,B,ans0
    end
    @inbounds while Index2 != Stop_criterion_bisect_step
        if Index1 == Box_dimension
            Index1 = 1
        else
            Index1 += 1
        end
        Size_stop=1
        for i=Index1:Index1+Box_dimension-1
            if B[Index_list[i]][2]-B[Index_list[i]][1] > Stop_criterion_subbox_size[Index_list[i]]
                Size_stop=0
                Index2 += 1
                B1 = SameBox(B)
                B1[Index_list[i]][2] = (B[Index_list[i]][1] + B[Index_list[i]][2]) / 2\n")
        if MC_method=="Simple"
            write(Temp_bisect_search_file,"                ans1 = sumo(B1,NN)\n")
        elseif MC_method=="Antithetic"
            write(Temp_bisect_search_file,"                ans1=sumo_antithetic(B1,NN)\n")
        end
        write(Temp_bisect_search_file,"                B2 = SameBox(B)
                B2[Index_list[i]][1] = (B[Index_list[i]][1] + B[Index_list[i]][2]) / 2\n")
        if MC_method=="Simple"
            write(Temp_bisect_search_file,"                ans2 = sumo(B2,NN)\n")
        elseif MC_method=="Antithetic"
            write(Temp_bisect_search_file,"                ans2=sumo_antithetic(B2,NN)\n")
        end
        write(Temp_bisect_search_file,"                integral_number+=2
                if ans1 > ans2
                    B=B1
                    ans0=ans1
                else
                    B=B2
                    ans0=ans2
                end
                if ans0 < 1.05
                    return \"No multistationary sub-box\",integral_number,Index2,B,ans0
                elseif ans0 > 2.95
                    return \"Found multistationary sub-box\",integral_number,Index2,B,ans0
                end
                break
            end
        end
        if Size_stop==1
            return \"Termination condition on the size of sub-boxes reached.\",integral_number,Index2,B,ans0
        end
    end
    return \"Termination condition on the number of bisect steps reached.\",integral_number,Index2,B,ans0
end\n")
    end
    close(Temp_bisect_search_file)
end

function Do_pretask(Input_lines)
    Do_subtask1(Input_lines)
    Do_subtask2()
end
