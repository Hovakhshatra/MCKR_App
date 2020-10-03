# The Simmple and Antithetic MCKR for example number 1 with parameters equipped with Uniform distibution.
# Here only T1 and T2 are free. And k1 until k6 are fixed to the following values;
# k=[0.7329,100,73.29,50,100,5]

import Distributions: Uniform # We need the "Uniform" function from the "Distribution" package.

function samplor(a, b) # To generate a uniform random number between a and b.
    return rand(Uniform(a, b))
end

function chi(a, b, x) # The indicator function on the interval [a,b].
    if x < a
        return 0
    elseif x > b
        return 0
    else
        return 1
    end
end

@fastmath function h(t)
   	return 366450 * t^2 + 537142.41 * t
end
@fastmath function beta(t)
   	return 2518322.5 * t^2 + 63502.1205 * t + 26857.1205
end
@fastmath function g1(t, T2) # g_{T_2}(t)
   	return beta(t) * (T2 - t) / h(t)
end
@fastmath function detJf(t, T1, T2)
   	return 537142.41 * T1 - 63502.1205 * T2 + 7554967.5 * t^2 + 2 * t * (366450 * T1 - 2518322.5 * T2 + 63502.1205) + 26857.1205
end
@fastmath function J(t, T2)
   	return abs(detJf(t, g1(t, T2), T2))
end

@fastmath function IntegrandWithoutCOEFF(t, aT1, bT1, T2)
   	return J(t, T2) * chi(aT1, bT1, g1(t, T2)) / h(t)
end

# MCKR simple. Input:Box, sample size. Output: I-hat.
@fastmath function sumo(B,NN)
    TT1 = B[1]
    TT2 = B[2]
    Ians = 0
    st1 = time_ns()
    COEFF = TT2[2] / (TT1[2] - TT1[1])
    t = samplor(0, TT2[2])
    T2 = samplor(TT2[1], TT2[2])
    Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2)
    Ians = COEFF * Ians
    @inbounds for i = 2:NN
        t = samplor(0, TT2[2])
   		T2 = samplor(TT2[1], TT2[2])
   		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
        Ians += delto / i
    end
   	return Ians
end

# MCKR simple. Input:Box, sample size. Output: I-hat, e-hat.
@fastmath function sumo_with_S(B,NN)
    TT1 = B[1]
    TT2 = B[2]
    Ians = 0
    S = 0
    st1 = time_ns()
    COEFF = TT2[2] / (TT1[2] - TT1[1])
    t = samplor(0, TT2[2])
    T2 = samplor(TT2[1], TT2[2])
    Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2)
    Ians = COEFF * Ians
    @inbounds for i = 2:NN
        t = samplor(0, TT2[2])
   		T2 = samplor(TT2[1], TT2[2])
   		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
        Ians += delto / i
      	S += delto^2 * (i - 1) / i
    end
   	return Ians, S
end

# MCKR simple. Input:Box, antithetic size. Output: I-hat.
@fastmath function sumo_antithetic(B,NN)
    TT1 = B[1]
    TT2 = B[2]
    Ians = 0
    centero = [(0 + TT2[2]) / 2,(TT2[1] + TT2[2]) / 2]
    COEFF = TT2[2] / (TT1[2] - TT1[1])
    t = samplor(0, TT2[2])
    T2 = samplor(TT2[1], TT2[2])
    Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2)
    Ians = COEFF * Ians
    t = 2 * centero[1] - t
   	T2 = 2 * centero[2] - T2
   	delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
    Ians += delto / 2
    @inbounds for i = 1:floor(NN / 2)
        t = samplor(0, TT2[2])
   		T2 = samplor(TT2[1], TT2[2])
   		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
        Ians += delto / (2 * i + 1)
        t = 2 * centero[1] - t
       	T2 = 2 * centero[2] - T2
        delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
        Ians += delto / (2 * i + 2)
    end
   	return Ians
end

# MCKR simple. Input:Box, antithetic size. Output: I-hat, e-hat.
@fastmath function sumo_antithetic_with_S(B,NN)
    TT1 = B[1]
    TT2 = B[2]
    Ians = 0
    S = 0
    centero = [(0 + TT2[2]) / 2,(TT2[1] + TT2[2]) / 2]
    COEFF = TT2[2] / (TT1[2] - TT1[1])
    t = samplor(0, TT2[2])
    T2 = samplor(TT2[1], TT2[2])
    Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2)
    Ians = COEFF * Ians
    t = 2 * centero[1] - t
   	T2 = 2 * centero[2] - T2
   	delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
    Ians += delto / 2
   	S += delto^2 / 2
    @inbounds for i = 1:floor(NN / 2)
        t = samplor(0, TT2[2])
   		T2 = samplor(TT2[1], TT2[2])
   		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
        Ians += delto / (2 * i + 1)
      	S += (delto^2) * 2 * i / (2 * i + 1)
        t = 2 * centero[1] - t
       	T2 = 2 * centero[2] - T2
        delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
        Ians += delto / (2 * i + 2)
      	S += (delto^2) * (2 * i + 1) / (2 * i + 2)
    end
   	return Ians, S
end
