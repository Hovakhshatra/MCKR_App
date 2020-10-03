# This is the HQ (head quarter) which reads the Input.txt file and choose the appropriate julia file.

User_input_file_name=ARGS[1] # get the input file name from the command line the same time of calling this file.

Input_file=open(joinpath(pwd(),"$User_input_file_name"),"r")
Input_lines=readlines(Input_file)
close(Input_file)

Report_lines=["\u066d\u066D\u066D\u066D\u066D Report \u066D\u066D\u066D\u066D\u066D\n","\n","\u066D User input\n"]

# Writing the input given by the user at the beginning of the report.
for i=1:length(Input_lines)
    push!(Report_lines,Input_lines[i]*"\n")
end

push!(Report_lines,"\u066D End of input.\n")
push!(Report_lines,"\n")

# Detecting the task.
for i=1:length(Input_lines)
    line=Input_lines[i]
    if length(line)>6 && line[1:6]=="Task: "
        if line[7:end]=="MCKR parallel compare"
            global Responsible_Julia_file_name="MCKR_parallel_compare"
            deleteat!(Input_lines,i)
            break
        elseif line[7:end]=="MCKR integrate"
            global Responsible_Julia_file_name="MCKR_integrate"
            deleteat!(Input_lines,i)
            break
        elseif line[7:end]=="MCKR integrate parallelized"
            global Responsible_Julia_file_name="MCKR_integrate_paral"
            deleteat!(Input_lines,i)
            break
        elseif line[7:end]=="MCKR integrate bisect search"
            global Responsible_Julia_file_name="MCKR_bisect_search"
            deleteat!(Input_lines,i)
            break
        elseif line[7:end]=="MCKR integrate bisect search parallelized"
            global Responsible_Julia_file_name="MCKR_bisect_search_paral"
            deleteat!(Input_lines,i)
            break
        elseif line[7:end]=="MCKR integrate bisect search parallelized - include tracking"
            global Responsible_Julia_file_name="MCKR_bisect_search_paral_tracker"
            deleteat!(Input_lines,i)
            break
        elseif line[7:end]=="MCKR method comparison table parallelized"
            global Responsible_Julia_file_name="MCKR_method_comparison_table_paral"
            deleteat!(Input_lines,i)
            break
        elseif line[7:end]=="MCKR method comparison table parallelized - include tracking"
            global Responsible_Julia_file_name="MCKR_method_comparison_table_paral_tracker"
            deleteat!(Input_lines,i)
            break
        end
    end
end

# Detecting the requested path to save the report file.
Save_location=pwd() # Default path.
for i=1:length(Input_lines)
    line=Input_lines[i]
    if length(line)>15 && line[1:15]=="Save_location: "
        global Save_location=line[16:end]
        deleteat!(Input_lines,i)
        break
    end
end

# Detecting the requested name for the report file.
Save_name="Output.txt"
for i=1:length(Input_lines)
    line=Input_lines[i]
    if length(line)>11 && line[1:11]=="Save_name: "
        global Save_name=line[12:end]*".txt"
        deleteat!(Input_lines,i)
        break
    end
end

# Detecting if the user wants the cpu information being included.
cpu_info_request=0
for i=1:length(Input_lines)
    line=Input_lines[i]
    if length(line)==23 && line[1:23]=="Include cpu information"
        global cpu_info_request=1
        deleteat!(Input_lines,i)
        break
    end
end

# Detecting if the user wants the memory information being included.
memory_info_request=0
for i=1:length(Input_lines)
    line=Input_lines[i]
    if length(line)==26 && line[1:26]=="Include memory information"
        global memory_info_request=1
        global memory_info_total=Int(Sys.total_memory())/2^30
        global memory_info_free=Int(Sys.free_memory())/2^30
        deleteat!(Input_lines,i)
        break
    end
end

# Doing the task.
include(joinpath(pwd(),"Tasks",Responsible_Julia_file_name*"_pretask.jl"))
Do_pretask(Input_lines)
include(joinpath(pwd(),"Tasks",Responsible_Julia_file_name*".jl"))
Task_result=Do_task()

# Adding the resport from the responsible Julia file to the report.
for i=1:length(Task_result)
    push!(Report_lines,Task_result[i])
end

# Adding the cpu information if it is asked.
if cpu_info_request==1
    push!(Report_lines,"\n\n\u066d Information of cpus of the computing system.\n")
    push!(Report_lines,"This system has $(length(Sys.cpu_info())) cpus.\n")
    push!(Report_lines,"$(Sys.cpu_info())\n")
end

# Adding the memory information if it is asked.
if memory_info_request==1
    push!(Report_lines,"\n\u066d Information about the memory.\n")
    push!(Report_lines,"The total memory is $memory_info_total GB.\n")
    push!(Report_lines,"The free memory before compiling the task was $memory_info_free GB.\n")
end

# Now writing the report text file.
Output_file=open(joinpath(Save_location,Save_name),"w")
for i=1:length(Report_lines)
    write(Output_file,Report_lines[i])
end
close(Output_file)
