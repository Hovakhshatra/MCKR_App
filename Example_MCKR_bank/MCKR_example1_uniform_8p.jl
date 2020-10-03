# The Simmple and Antithetic MCKR for example number 1 with parameters equipped with Uniform distibution.
# All 8 parameters are free.

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

@fastmath function h1(t, k1, k2, k3, k4, k5, k6)
   	return(k1 * k2 * k4 * k5 * t^2 + k1 * k2 * k3 * k5 * t)
end
@fastmath function beta1(t, k1, k2, k3, k4, k5, k6)
   	return((k4 * k5 * (k1 + k2) * t^2 + k1 * k5 * (k2 + k3) * t + k1 * k2 * k3) * k6)
end
@fastmath function g1(t, T2, k1, k2, k3, k4, k5, k6)
   	return beta1(t, k1, k2, k3, k4, k5, k6) * (T2 - t) / h1(t, k1, k2, k3, k4, k5, k6)
end
@fastmath function detJf(t, T1, T2, k1, k2, k3, k4, k5, k6)
   	return(T1 * k1 * k2 * k3 * k5 - T2 * k1 * k2 * k5 * k6 - T2 * k1 * k3 * k5 * k6 + k1 * k2 * k3 * k6 + 3 * t^2 * (k1 * k4 * k5 * k6
	       + k2 * k4 * k5 * k6) + 2 * t * (T1 * k1 * k2 * k4 * k5 - T2 * k1 * k4 * k5 * k6 - T2 * k2 * k4 * k5 * k6 + k1 * k2 * k5 * k6 + k1 * k3 * k5 * k6))
end
@fastmath function J(t, T2, k1, k2, k3, k4, k5, k6)
   	return abs(detJf(t, g1(t, T2, k1, k2, k3, k4, k5, k6), T2, k1, k2, k3, k4, k5, k6))
end

@fastmath function IntegrandWithoutCOEFF(t, aT1, bT1, T2, k1, k2, k3, k4, k5, k6)
   	return J(t, T2, k1, k2, k3, k4, k5, k6) * chi(aT1, bT1, g1(t, T2, k1, k2, k3, k4, k5, k6)) / h1(t, k1, k2, k3, k4, k5, k6)
end

# MCKR simple. Input:Box, sample size. Output: I-hat.
@fastmath function sumo(B, NN)
	kk1=B[1]
	kk2=B[2]
	kk3=B[3]
	kk4=B[4]
	kk5=B[5]
	kk6=B[6]
	TT1=B[7]
	TT2=B[8]
   	Ians = 0
   	COEFF = TT2[2] / (TT1[2] - TT1[1])
   	t = samplor(0, TT2[2])
   	k1 = samplor(kk1[1], kk1[2])
   	k2 = samplor(kk2[1], kk2[2])
   	k3 = samplor(kk3[1], kk3[2])
   	k4 = samplor(kk4[1], kk4[2])
   	k5 = samplor(kk5[1], kk5[2])
   	k6 = samplor(kk6[1], kk6[2])
   	T2 = samplor(TT2[1], TT2[2])
   	Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6)
   	Ians = COEFF * Ians
   	@inbounds for i = 2:NN
      		t = samplor(0, TT2[2])
      		k1 = samplor(kk1[1], kk1[2])
      		k2 = samplor(kk2[1], kk2[2])
      		k3 = samplor(kk3[1], kk3[2])
      		k4 = samplor(kk4[1], kk4[2])
      		k5 = samplor(kk5[1], kk5[2])
      		k6 = samplor(kk6[1], kk6[2])
      		T2 = samplor(TT2[1], TT2[2])
      		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
      		Ians += delto / i
   	end
   	return Ians
end

# MCKR simple. Input:Box, sample size. Output: I-hat, e-hat.
@fastmath function sumo_with_S(B, NN)
	kk1=B[1]
	kk2=B[2]
	kk3=B[3]
	kk4=B[4]
	kk5=B[5]
	kk6=B[6]
	TT1=B[7]
	TT2=B[8]
   	S = 0
   	Ians = 0
   	COEFF = TT2[2] / (TT1[2] - TT1[1])
   	t = samplor(0, TT2[2])
   	k1 = samplor(kk1[1], kk1[2])
   	k2 = samplor(kk2[1], kk2[2])
   	k3 = samplor(kk3[1], kk3[2])
   	k4 = samplor(kk4[1], kk4[2])
   	k5 = samplor(kk5[1], kk5[2])
   	k6 = samplor(kk6[1], kk6[2])
   	T2 = samplor(TT2[1], TT2[2])
   	Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6)
   	Ians = COEFF * Ians
   	@inbounds for i = 2:NN
      		t = samplor(0, TT2[2])
      		k1 = samplor(kk1[1], kk1[2])
      		k2 = samplor(kk2[1], kk2[2])
      		k3 = samplor(kk3[1], kk3[2])
      		k4 = samplor(kk4[1], kk4[2])
      		k5 = samplor(kk5[1], kk5[2])
      		k6 = samplor(kk6[1], kk6[2])
      		T2 = samplor(TT2[1], TT2[2])
      		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
      		Ians += delto / i
      		S += delto^2 * (i - 1) / i
   	end
   	return Ians, S
end

# MCKR antithetic. Input:Box, sample size. Output: I-hat.
@fastmath function sumo_antithetic(B, NN)
	kk1=B[1]
	kk2=B[2]
	kk3=B[3]
	kk4=B[4]
	kk5=B[5]
	kk6=B[6]
	TT1=B[7]
	TT2=B[8]
   	Ians = 0
   	centero = [(0 + TT2[2]) / 2,(kk1[1] + kk1[2]) / 2,(kk2[1] + kk2[2]) / 2,(kk3[1] + kk3[2]) / 2,(kk4[1] + kk4[2]) / 2,(kk5[1] + kk5[2]) / 2,(kk6[1] + kk6[2]) / 2
	        ,(TT2[1] + TT2[2]) / 2]
   	COEFF = TT2[2] / (TT1[2] - TT1[1])
   	t = samplor(0, TT2[2])
   	k1 = samplor(kk1[1], kk1[2])
   	k2 = samplor(kk2[1], kk2[2])
   	k3 = samplor(kk3[1], kk3[2])
   	k4 = samplor(kk4[1], kk4[2])
   	k5 = samplor(kk5[1], kk5[2])
   	k6 = samplor(kk6[1], kk6[2])
   	T2 = samplor(TT2[1], TT2[2])
   	Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6)
   	Ians = COEFF * Ians
   	t = 2 * centero[1] - t
   	k1 = 2 * centero[2] - k1
   	k2 = 2 * centero[3] - k2
   	k3 = 2 * centero[4] - k3
   	k4 = 2 * centero[5] - k4
   	k5 = 2 * centero[6] - k5
   	k6 = 2 * centero[7] - k6
   	T2 = 2 * centero[8] - T2
   	delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
   	Ians += delto / 2
   	@inbounds for i = 1:floor(NN / 2)
      		t = samplor(0, TT2[2])
      		k1 = samplor(kk1[1], kk1[2])
      		k2 = samplor(kk2[1], kk2[2])
      		k3 = samplor(kk3[1], kk3[2])
      		k4 = samplor(kk4[1], kk4[2])
      		k5 = samplor(kk5[1], kk5[2])
      		k6 = samplor(kk6[1], kk6[2])
      		T2 = samplor(TT2[1], TT2[2])
      		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
      		Ians += delto / (2 * i + 1)
      		t = 2 * centero[1] - t
      		k1 = 2 * centero[2] - k1
      		k2 = 2 * centero[3] - k2
      		k3 = 2 * centero[4] - k3
      		k4 = 2 * centero[5] - k4
      		k5 = 2 * centero[6] - k5
      		k6 = 2 * centero[7] - k6
      		T2 = 2 * centero[8] - T2
      		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
      		Ians += delto / (2 * i + 2)
   	end
   	return Ians
end

# MCKR antithetic. Input:Box, sample size. Output: I-hat, e-hat.
@fastmath function sumo_antithetic_with_S(B, NN)
	kk1=B[1]
	kk2=B[2]
	kk3=B[3]
	kk4=B[4]
	kk5=B[5]
	kk6=B[6]
	TT1=B[7]
	TT2=B[8]
   	S = 0
   	Ians = 0
   	centero = [(0 + TT2[2]) / 2,(kk1[1] + kk1[2]) / 2,(kk2[1] + kk2[2]) / 2,(kk3[1] + kk3[2]) / 2,(kk4[1] + kk4[2]) / 2,(kk5[1] + kk5[2]) / 2,(kk6[1] + kk6[2]) / 2
	        ,(TT2[1] + TT2[2]) / 2]
   	COEFF = TT2[2] / (TT1[2] - TT1[1])
   	t = samplor(0, TT2[2])
   	k1 = samplor(kk1[1], kk1[2])
   	k2 = samplor(kk2[1], kk2[2])
   	k3 = samplor(kk3[1], kk3[2])
   	k4 = samplor(kk4[1], kk4[2])
   	k5 = samplor(kk5[1], kk5[2])
   	k6 = samplor(kk6[1], kk6[2])
   	T2 = samplor(TT2[1], TT2[2])
   	Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6)
   	Ians = COEFF * Ians
   	t = 2 * centero[1] - t
   	k1 = 2 * centero[2] - k1
   	k2 = 2 * centero[3] - k2
   	k3 = 2 * centero[4] - k3
   	k4 = 2 * centero[5] - k4
   	k5 = 2 * centero[6] - k5
   	k6 = 2 * centero[7] - k6
   	T2 = 2 * centero[8] - T2
   	delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
   	Ians += delto / 2
   	S += delto^2 / 2
   	@inbounds for i = 1:floor(NN / 2)
      		t = samplor(0, TT2[2])
      		k1 = samplor(kk1[1], kk1[2])
      		k2 = samplor(kk2[1], kk2[2])
      		k3 = samplor(kk3[1], kk3[2])
      		k4 = samplor(kk4[1], kk4[2])
      		k5 = samplor(kk5[1], kk5[2])
      		k6 = samplor(kk6[1], kk6[2])
      		T2 = samplor(TT2[1], TT2[2])
      		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
      		Ians += delto / (2 * i + 1)
      		S += (delto^2) * 2 * i / (2 * i + 1)
      		t = 2 * centero[1] - t
      		k1 = 2 * centero[2] - k1
      		k2 = 2 * centero[3] - k2
      		k3 = 2 * centero[4] - k3
      		k4 = 2 * centero[5] - k4
      		k5 = 2 * centero[6] - k5
      		k6 = 2 * centero[7] - k6
      		T2 = 2 * centero[8] - T2
      		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
      		Ians += delto / (2 * i + 2)
      		S += (delto^2) * (2 * i + 1) / (2 * i + 2)
   	end
   	return Ians, S
end
