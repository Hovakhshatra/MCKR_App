# The Simmple and Antithetic MCKR for example number 2 with parameters equipped with normal distibution.

# We need to define two functions TNpdf and samplorTN in below.
using Distributions
function samplor(a, b)
   	return rand(Uniform(a, b))
end
function TNpdf(a, b, m, v, x)
   	return(Distributions.pdf(Distributions.Truncated(Normal(m, v), a, b), x))
end
function samplorTN(a, b, m, v)
   	return(rand(Distributions.Truncated(Normal(m, v), a, b)))
end

@fastmath function alok1(x1,x2,k2,k3,k4)
	return (2*k2*x1^3+k3*x1*x2^2-2*k4*x2^3)/(x2*x1^2)
end
@fastmath function aloT(x1,x2,k2,k3,k4)
	return x1+x2
end
@fastmath function detJf(x1,x2,k1,k2,k3,k4)
	return -(6*k2+k1)*x1^2-(6*k4+k3)*x2^2+2*(k1+k3)*x1*x2
end
@fastmath function J(x1,x2,k2,k3,k4)
	return abs(detJf(x1,x2,alok1(x1,x2,k2,k3,k4),k2,k3,k4))
end

@fastmath function IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, ak1, bk1, mk1, vk1, aT, bT, mT, vT)
   	return J(x1, x2, k2, k3, k4) * TNpdf(ak1, bk1, mk1, vk1, alok1(x1, x2, k2, k3, k4)) * TNpdf(aT, bT, mT, vT, aloT(x1, x2, k2, k3, k4)) / (x2 * x1^2)
end

# MCKR simple. Input:Box, sample size. Output: I-hat.
@fastmath function sumo(B, NN) # Each kk here is a 4-tuple; a_i, b_i, mu_i, sigma_i^2.
	kk1=B[1]
	kk2=B[2]
	kk3=B[3]
	kk4=B[4]
	TT=B[5]
   	Ians = 0
   	COEFF = (TT[2]^2)
   	x1 = samplor(0, TT[2])
   	x2 = samplor(0, TT[2])
   	k2 = samplorTN(kk2[1], kk2[2], kk2[3], kk2[4])
   	k3 = samplorTN(kk3[1], kk3[2], kk3[3], kk3[4])
   	k4 = samplorTN(kk4[1], kk4[2], kk4[3], kk4[4])
   	Ians += IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4])
   	Ians = COEFF * Ians
   	@inbounds for i = 2:NN
      		x1 = samplor(0, TT[2])
      		x2 = samplor(0, TT[2])
      		k2 = samplorTN(kk2[1], kk2[2], kk2[3], kk2[4])
      		k3 = samplorTN(kk3[1], kk3[2], kk3[3], kk3[4])
      		k4 = samplorTN(kk4[1], kk4[2], kk4[3], kk4[4])
      		delto = COEFF * IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4]) - Ians
      		Ians += delto / (i + 1)
   	end
   	return Ians
end

# MCKR simple. Input:Box, sample size. Output: I-hat, e-hat.
@fastmath function sumo_with_S(B, NN) # Each kk here is a 4-tuple; a_i, b_i, mu_i, sigma_i^2.
	kk1=B[1]
	kk2=B[2]
	kk3=B[3]
	kk4=B[4]
	TT=B[5]
	S = 0
   	Ians = 0
   	COEFF = (TT[2]^2)
   	x1 = samplor(0, TT[2])
   	x2 = samplor(0, TT[2])
   	k2 = samplorTN(kk2[1], kk2[2], kk2[3], kk2[4])
   	k3 = samplorTN(kk3[1], kk3[2], kk3[3], kk3[4])
   	k4 = samplorTN(kk4[1], kk4[2], kk4[3], kk4[4])
   	Ians += IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4])
   	Ians = COEFF * Ians
   	@inbounds for i = 2:NN
      		x1 = samplor(0, TT[2])
      		x2 = samplor(0, TT[2])
      		k2 = samplorTN(kk2[1], kk2[2], kk2[3], kk2[4])
      		k3 = samplorTN(kk3[1], kk3[2], kk3[3], kk3[4])
      		k4 = samplorTN(kk4[1], kk4[2], kk4[3], kk4[4])
      		delto = COEFF * IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4]) - Ians
      		Ians += delto / (i + 1)
      		S += (delto^2) * i / (i + 1)
   	end
   	return Ians, S
end

# MCKR antithetic. Input:Box, antithetic size. Output: I-hat.
@fastmath function sumo_antithetic(B, NN)
	kk1=B[1]
	kk2=B[2]
	kk3=B[3]
	kk4=B[4]
	TT=B[5]
   	Ians = 0
   	if (kk2[1] + kk2[2]) / 2 != kk2[3] || (kk3[1] + kk3[2]) / 2 != kk3[3] || (kk4[1] + kk4[2]) / 2 != kk4[3]
      		return error # Because for antithetic Monte-Carlo we need symmetric distribtuion (with respect to its support, not only the function alone).
   	end
   	centero = [(0 + TT[2]) / 2,(0 + TT[2]) / 2,kk2[3],kk3[3],kk4[3]]
   	COEFF = (TT[2]^2)
   	x1 = samplor(0, TT[2])
   	x2 = samplor(0, TT[2])
   	k2 = samplorTN(kk2[1], kk2[2], kk2[3], kk2[4])
   	k3 = samplorTN(kk3[1], kk3[2], kk3[3], kk3[4])
   	k4 = samplorTN(kk4[1], kk4[2], kk4[3], kk4[4])
   	Ians += IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4])
   	Ians = COEFF * Ians
   	x1 = 2 * centero[1] - x1
   	x2 = 2 * centero[2] - x2
   	k2 = 2 * centero[3] - k2
   	k3 = 2 * centero[4] - k3
   	k4 = 2 * centero[5] - k4
   	delto = COEFF * IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4]) - Ians
   	Ians += delto / 2
   	@inbounds for i = 1:floor(NN / 2)
      		x1 = samplor(0, TT[2])
      		x2 = samplor(0, TT[2])
      		k2 = samplorTN(kk2[1], kk2[2], kk2[3], kk2[4])
      		k3 = samplorTN(kk3[1], kk3[2], kk3[3], kk3[4])
      		k4 = samplorTN(kk4[1], kk4[2], kk4[3], kk4[4])
      		delto = COEFF * IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4]) - Ians
      		Ians += delto / (2 * i + 1)
      		x1 = 2 * centero[1] - x1
      		x2 = 2 * centero[2] - x2
      		k2 = 2 * centero[3] - k2
      		k3 = 2 * centero[4] - k3
      		k4 = 2 * centero[5] - k4
      		delto = COEFF * IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4]) - Ians
      		Ians += delto / (2 * i + 2)
   	end
   	return Ians
end

# MCKR antithetic. Input:Box, antithetic size. Output: I-hat, e-hat.
@fastmath function sumo_antithetic_with_S(B, NN)
	kk1=B[1]
	kk2=B[2]
	kk3=B[3]
	kk4=B[4]
	TT=B[5]
	S = 0
   	Ians = 0
   	if (kk2[1] + kk2[2]) / 2 != kk2[3] || (kk3[1] + kk3[2]) / 2 != kk3[3] || (kk4[1] + kk4[2]) / 2 != kk4[3]
      		return error # Because for antithetic Monte-Carlo we need symmetric distribtuion (with respect to its support, not only the function alone).
   	end
   	centero = [(0 + TT[2]) / 2,(0 + TT[2]) / 2,kk2[3],kk3[3],kk4[3]]
   	COEFF = (TT[2]^2)
   	x1 = samplor(0, TT[2])
   	x2 = samplor(0, TT[2])
   	k2 = samplorTN(kk2[1], kk2[2], kk2[3], kk2[4])
   	k3 = samplorTN(kk3[1], kk3[2], kk3[3], kk3[4])
   	k4 = samplorTN(kk4[1], kk4[2], kk4[3], kk4[4])
   	Ians += IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4])
   	Ians = COEFF * Ians
   	x1 = 2 * centero[1] - x1
   	x2 = 2 * centero[2] - x2
   	k2 = 2 * centero[3] - k2
   	k3 = 2 * centero[4] - k3
   	k4 = 2 * centero[5] - k4
   	delto = COEFF * IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4]) - Ians
   	Ians += delto / 2
   	S += (delto^2) * 1 / 2
   	@inbounds for i = 1:floor(NN / 2)
      		x1 = samplor(0, TT[2])
      		x2 = samplor(0, TT[2])
      		k2 = samplorTN(kk2[1], kk2[2], kk2[3], kk2[4])
      		k3 = samplorTN(kk3[1], kk3[2], kk3[3], kk3[4])
      		k4 = samplorTN(kk4[1], kk4[2], kk4[3], kk4[4])
      		delto = COEFF * IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4]) - Ians
      		Ians += delto / (2 * i + 1)
      		S += (delto^2) * 2 * i / (2 * i + 1)
      		x1 = 2 * centero[1] - x1
      		x2 = 2 * centero[2] - x2
      		k2 = 2 * centero[3] - k2
      		k3 = 2 * centero[4] - k3
      		k4 = 2 * centero[5] - k4
      		delto = COEFF * IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4]) - Ians
      		Ians += delto / (2 * i + 2)
      		S += (delto^2) * (2 * i + 1) / (2 * i + 2)
   	end
   	return Ians,S
end
