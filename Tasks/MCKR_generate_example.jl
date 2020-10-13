# The MCKR generate example task.
# This MCKR task is going to generate a new file in the example bank folder
# corresponding to the parametric system of equations given by the user.

include(joinpath(pwd(),"Temp","MCKR_generate_example_input.jl"))

function Do_task()
    # The following string is going to be used later to check symmetry of the distribution for the normal case.
    Disjunction_string=""
    for i=1:Number_of_parameters
        Disjunction_string=Disjunction_string*"(B[$i][1]+B[$i][2])/2!=B[$i][3] || "
    end
    Disjunction_string="if "*Disjunction_string[1:end-4]*"\n    return error # Because for antithetic Monte-Carlo we need symmetric distribtuion (with respect to its support, not only the function alone).
end\n    "
    # Pieces of codes different for Uniform and Normal distributions.
    Distributions_dict=Dict(
        "Uniform"=>["import Distributions: Uniform
function samplor(a, b)
    return rand(Uniform(a, b))
end
function chi(a, b, x)
    if x < a
        return 0
    elseif x > b
        return 0
    else
        return 1
    end
end","samplor",""],
        "Normal"=>["using Distributions
function samplor(a, b)
   	return rand(Uniform(a, b))
end
function TNpdf(a, b, m, v, x)
   	return(Distributions.pdf(Distributions.Truncated(Normal(m, v), a, b), x))
end
function samplorTN(a, b, m, v)
   	return(rand(Distributions.Truncated(Normal(m, v), a, b)))
end","samplorTN",Disjunction_string]
    )
    # Some strings that get repeated very often.
    arguments_string=""
    for j=1:Number_of_variables
        arguments_string=arguments_string*"x$j,"
    end
    for j=Number_of_variables+1:Number_of_parameters-1
        arguments_string=arguments_string*"k$j,"
    end
    arguments_string=arguments_string*"k$Number_of_parameters"
    # Creating the example file.
    Example_file=open(joinpath(pwd(),"Example_MCKR_bank","$Example_name"),"w")
    write(Example_file,"$(Distributions_dict[Distribution][1])\n")
    for i=1:Number_of_variables
        write(Example_file,"@fastmath function g$i($arguments_string)
    return $(g_dict[i])
end\n")
    end
    write(Example_file,"@fastmath function J($arguments_string)
    return abs($Determinant_expression)
end
@fastmath function s($arguments_string,B)
    return J($arguments_string)*")
    if Distribution=="Uniform"
        for i=1:Number_of_variables-1
            write(Example_file,"chi(B[$i][1],B[$i][2],g$i($arguments_string))*")
        end
        write(Example_file,"chi(B[$Number_of_variables][1],B[$Number_of_variables][2],g$Number_of_variables($arguments_string))\n")
    elseif Distribution=="Normal"
        for i=1:Number_of_variables-1
            write(Example_file,"TNpdf(B[$i][1],B[$i][2],B[$i][3],B[$i][4],g$i($arguments_string))*")
        end
        write(Example_file,"TNpdf(B[$Number_of_variables][1],B[$Number_of_variables][2],B[$Number_of_variables][3],B[$Number_of_variables][4],g$Number_of_variables($arguments_string))\n")
    end
    write(Example_file,"end
@fastmath function IntegrandWithoutCOEFF($arguments_string,B)
    return $summand_expression
end

# MCKR simple. Input:Box, sample size. Output: I-hat.
@fastmath function sumo(B,NN)
    Ians=0
    COEFF=")
    COEFF_string="1"
    if Distribution=="Uniform"
        COEFF_string=COEFF_string*"/("
        for i=1:Number_of_variables
            COEFF_string=COEFF_string*"(B[$i][2]-B[$i][1])*"
        end
        COEFF_string=COEFF_string[1:end-1]*")"
    end
    write(Example_file,"$COEFF_string\n")
    for i=1:Number_of_variables
        write(Example_file,"    x$i=samplor(0,1)\n")
    end
    for i=Number_of_variables+1:Number_of_parameters
        if Distribution=="Uniform"
            write(Example_file,"    k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2])\n")
        elseif Distribution=="Normal"
            write(Example_file,"    k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2],B[$i][3],B[$i][4])\n")
        end
    end
    write(Example_file,"    Ians+=IntegrandWithoutCOEFF($arguments_string,B)
    Ians=COEFF*Ians
    @inbounds for i=2:NN\n")
    for i=1:Number_of_variables
        write(Example_file,"        x$i=samplor(0,1)\n")
    end
    for i=Number_of_variables+1:Number_of_parameters
        if Distribution=="Uniform"
            write(Example_file,"        k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2])\n")
        elseif Distribution=="Normal"
            write(Example_file,"        k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2],B[$i][3],B[$i][4])\n")
        end
    end
    write(Example_file,"        delto=COEFF*IntegrandWithoutCOEFF($arguments_string,B)-Ians
        Ians+=delto/(i+1)
    end
    return Ians
end

# MCKR simple. Input:Box, sample size. Output: I-hat, e-hat.
@fastmath function sumo_with_S(B,NN)
    S=0
    Ians=0
    COEFF=$COEFF_string\n")
    for i=1:Number_of_variables
        write(Example_file,"    x$i=samplor(0,1)\n")
    end
    for i=Number_of_variables+1:Number_of_parameters
        if Distribution=="Uniform"
            write(Example_file,"    k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2])\n")
        elseif Distribution=="Normal"
            write(Example_file,"    k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2],B[$i][3],B[$i][4])\n")
        end
    end
    write(Example_file,"    Ians+=IntegrandWithoutCOEFF($arguments_string,B)
    Ians=COEFF*Ians
    @inbounds for i=2:NN\n")
    for i=1:Number_of_variables
        write(Example_file,"        x$i=samplor(0,1)\n")
    end
    for i=Number_of_variables+1:Number_of_parameters
        if Distribution=="Uniform"
            write(Example_file,"        k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2])\n")
        elseif Distribution=="Normal"
            write(Example_file,"        k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2],B[$i][3],B[$i][4])\n")
        end
    end
    write(Example_file,"        delto=COEFF*IntegrandWithoutCOEFF($arguments_string,B)-Ians
        Ians+=delto/(i+1)
        S+=(delto^2)*i/(i+1)
    end
    return Ians,S
end

# MCKR antithetic. Input:Box, antithetic size. Output: I-hat.
@fastmath function sumo_antithetic(B,NN)
    Ians=0
    $(Distributions_dict[Distribution][3])centero=")
    centero_string="["
    for i=1:Number_of_variables
        centero_string=centero_string*"(0+1)/2,"
    end
    for i=Number_of_variables+1:Number_of_parameters
        centero_string=centero_string*"(B[$i][1]+B[$i][2])/2,"
    end
    centero_string=centero_string[1:end-1]*"]"
    write(Example_file,"$centero_string
    COEFF=$COEFF_string\n")
    for i=1:Number_of_variables
        write(Example_file,"    x$i=samplor(0,1)\n")
    end
    for i=Number_of_variables+1:Number_of_parameters
        if Distribution=="Uniform"
            write(Example_file,"    k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2])\n")
        elseif Distribution=="Normal"
            write(Example_file,"    k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2],B[$i][3],B[$i][4])\n")
        end
    end
    write(Example_file,"    Ians+=IntegrandWithoutCOEFF($arguments_string,B)
    Ians=COEFF*Ians\n")
    for i=1:Number_of_variables
        write(Example_file,"    x$i=2*centero[$i]-x$i\n")
    end
    for i=Number_of_variables+1:Number_of_parameters
        write(Example_file,"    k$i=2*centero[$i]-k$i\n")
    end
    write(Example_file,"    delto=COEFF*IntegrandWithoutCOEFF($arguments_string,B)-Ians
    Ians+=delto/2
    @inbounds for i=1:floor(NN/2)\n")
    for i=1:Number_of_variables
        write(Example_file,"        x$i=samplor(0,1)\n")
    end
    for i=Number_of_variables+1:Number_of_parameters
        if Distribution=="Uniform"
            write(Example_file,"        k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2])\n")
        elseif Distribution=="Normal"
            write(Example_file,"        k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2],B[$i][3],B[$i][4])\n")
        end
    end
    write(Example_file,"        delto=COEFF*IntegrandWithoutCOEFF($arguments_string,B)-Ians
        Ians+=delto/(2*i+1)\n")
    for i=1:Number_of_variables
        write(Example_file,"        x$i=2*centero[$i]-x$i\n")
    end
    for i=Number_of_variables+1:Number_of_parameters
        write(Example_file,"        k$i=2*centero[$i]-k$i\n")
    end
    write(Example_file,"        delto=COEFF*IntegrandWithoutCOEFF($arguments_string,B)-Ians
        Ians+=delto/(2*i+2)
    end
    return Ians
end

# MCKR antithetic. Input:Box, antithetic size. Output: I-hat, e-hat.
@fastmath function sumo_antithetic_with_S(B,NN)
    S=0
    Ians=0
    $(Distributions_dict[Distribution][3])centero=")
    centero_string="["
    for i=1:Number_of_variables
        centero_string=centero_string*"(0+1)/2,"
    end
    for i=Number_of_variables+1:Number_of_parameters
        centero_string=centero_string*"(B[$i][1]+B[$i][2])/2,"
    end
    centero_string=centero_string[1:end-1]*"]"
    write(Example_file,"$centero_string
    COEFF=$COEFF_string\n")
    for i=1:Number_of_variables
        write(Example_file,"    x$i=samplor(0,1)\n")
    end
    for i=Number_of_variables+1:Number_of_parameters
        if Distribution=="Uniform"
            write(Example_file,"    k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2])\n")
        elseif Distribution=="Normal"
            write(Example_file,"    k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2],B[$i][3],B[$i][4])\n")
        end
    end
    write(Example_file,"    Ians+=IntegrandWithoutCOEFF($arguments_string,B)
    Ians=COEFF*Ians\n")
    for i=1:Number_of_variables
        write(Example_file,"    x$i=2*centero[$i]-x$i\n")
    end
    for i=Number_of_variables+1:Number_of_parameters
        write(Example_file,"    k$i=2*centero[$i]-k$i\n")
    end
    write(Example_file,"    delto=COEFF*IntegrandWithoutCOEFF($arguments_string,B)-Ians
    Ians+=delto/2
   	S+=(delto^2)*1/2
    @inbounds for i=1:floor(NN/2)\n")
    for i=1:Number_of_variables
        write(Example_file,"        x$i=samplor(0,1)\n")
    end
    for i=Number_of_variables+1:Number_of_parameters
        if Distribution=="Uniform"
            write(Example_file,"        k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2])\n")
        elseif Distribution=="Normal"
            write(Example_file,"        k$i=$(Distributions_dict[Distribution][2])(B[$i][1],B[$i][2],B[$i][3],B[$i][4])\n")
        end
    end
    write(Example_file,"        delto=COEFF*IntegrandWithoutCOEFF($arguments_string,B)-Ians
        Ians+=delto/(2*i+1)
        S+=(delto^2)*2*i/(2*i+1)\n")
    for i=1:Number_of_variables
        write(Example_file,"        x$i=2*centero[$i]-x$i\n")
    end
    for i=Number_of_variables+1:Number_of_parameters
        write(Example_file,"        k$i=2*centero[$i]-k$i\n")
    end
    write(Example_file,"        delto=COEFF*IntegrandWithoutCOEFF($arguments_string,B)-Ians
        Ians+=delto/(2*i+2)
        S +=(delto^2)*(2*i+1)/(2*i+2)
    end
    return Ians,S
end")
    close(Example_file)
    # Returning the success message.
    print("Example file $Example_name is generated. You can find it in the Example_MCKR_bank folder.\n")
    return(["Example file $Example_name is generated. You can find it in the Example_MCKR_bank folder.\n"])
end
