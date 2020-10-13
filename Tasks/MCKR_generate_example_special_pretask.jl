# The MCKR generate example special pretask.
# This MCKR task is going to generate a new file in the example bank folder
# corresponding to the parametric system of equations given by the user.

function Do_subtask1(Input_lines)
    # Reading the inputs from the Input_lines and writes a temporary input julia file for the current tasks.
    Temp_input_file=open(joinpath(pwd(),"Temp","MCKR_generate_example_input.jl"),"w")
    Example_name="User_made_example.jl"
    for i=1:length(Input_lines) # Detecting the given name for the example by the user.
        line=Input_lines[i]
        if length(line)>14 && line[1:14]=="Example_name: "
            write(Temp_input_file,"Example_name=\"$(line[15:end]).jl\"\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    for i=1:length(Input_lines) # Detecting number of variables.
        line=Input_lines[i]
        if length(line)>21 && line[1:21]=="Number_of_variables: "
            write(Temp_input_file,"Number_of_variables=$(line[22:end])\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    for i=1:length(Input_lines) # Detecting number of parameters.
        line=Input_lines[i]
        if length(line)>22 && line[1:22]=="Number_of_parameters: "
            write(Temp_input_file,"Number_of_parameters=$(line[23:end])\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    for i=1:length(Input_lines) # Detecting the determinant expression.
        line=Input_lines[i]
        if length(line)>24 && line[1:24]=="Determinant_expression: "
            write(Temp_input_file,"Determinant_expression=\"$(line[25:end])\"\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    # Making a dict containing the g expressions.
    g_dict_string="g_dict=Dict("
    lines_to_delete=[]
    g_index=1
    for i=1:length(Input_lines) # Detect a new g expression if exists any.
        line=Input_lines[i]
        if length(line)>14 && line[1:14]=="g_expression: "
            g_dict_string=g_dict_string*"\n    $g_index=>\"$(line[15:end])\","
            append!(lines_to_delete,i)
            g_index+=1
        end
    end
    g_dict_string=g_dict_string[1:end-1]*"\n)\n"
    write(Temp_input_file,g_dict_string)
    for i=1:length(lines_to_delete)
        deleteat!(Input_lines,i)
    end
    Distribution="Uniform"
    for i=1:length(Input_lines) # Detecting the distributions on the parameters.
        line=Input_lines[i]
        if length(line)>14 && line[1:14]=="Distribution: "
            write(Temp_input_file,"Distribution=\"$(line[15:end])\"\n")
            deleteat!(Input_lines,i)
            break
        end
    end
    for i=1:length(Input_lines) # Detect the upper bound input for the variables.
        line=Input_lines[i]
        if length(line)>23 && line[1:23]=="Variables_upper_bound: "
            write(Temp_input_file,"Variables_upper_bound=[$(line[24:end])]")
            append!(lines_to_delete,i)
            g_index+=1
        end
    end
    close(Temp_input_file)
end

function Do_subtask2()
    include(joinpath(pwd(),"Temp","MCKR_generate_example_input.jl"))
    Temp_input_file=open(joinpath(pwd(),"Temp","MCKR_variables_upper_bound_dict.jl"),"w")
    Variables_upper_bound_string="Variables_upper_bound_dict=Dict("
    for i=1:Number_of_variables
        if Variables_upper_bound[i]==0
            Variables_upper_bound_string=Variables_upper_bound_string*"$i=>\"1\","
        else
            Variables_upper_bound_string=Variables_upper_bound_string*"$i=>\"B[$(Variables_upper_bound[i])][2]\","
        end
    end
    Variables_upper_bound_string=Variables_upper_bound_string[1:end-1]*")"
    write(Temp_input_file,Variables_upper_bound_string)
    close(Temp_input_file)
end

function Do_pretask(Input_lines)
    Do_subtask1(Input_lines)
    Do_subtask2()
end
