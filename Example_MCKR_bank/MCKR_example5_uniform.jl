# The Simmple and Antithetic MCKR for example number 5 with parameters equipped with Uniform distibution.

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

@fastmath function g(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k13,k14)
    return -(t^5*k1*k4*k5*k6*k7*k10*k11 + t^5*k1*k4*k5*k6*k8*k10*k11 + t^5*k2*k4*k5*k6*k7*k10*k11 + t^5*k2*k4*k5*k6*k8*k10*k11 - t^4*k1*k4*k5*k6*k7*k10*k11*k13 - t^4*k1*k4*k5*k6*k8*k10*k11*k13 + t^4*k1*k4*k5*k7*k8*k10*k11*k14 - t^4*k2*k4*k5*k6*k7*k10*k11*k13 - t^4*k2*k4*k5*k6*k8*k10*k11*k13 + t^4*k2*k4*k5*k7*k8*k10*k11*k14 + t^4*k1*k2*k5*k6*k7*k10*k11 + t^4*k1*k2*k5*k6*k8*k10*k11 + t^4*k1*k3*k5*k6*k7*k10*k11 + t^4*k1*k3*k5*k6*k8*k10*k11 + t^4*k1*k4*k5*k6*k7*k8*k11 + t^4*k1*k4*k5*k6*k7*k9*k11 + t^4*k2*k4*k5*k6*k7*k8*k11 + t^4*k2*k4*k5*k6*k7*k9*k11 - t^3*k1*k2*k5*k6*k7*k10*k11*k13 - t^3*k1*k2*k5*k6*k8*k10*k11*k13 + t^3*k1*k2*k5*k7*k8*k10*k11*k14 - t^3*k1*k3*k5*k6*k7*k10*k11*k13 - t^3*k1*k3*k5*k6*k8*k10*k11*k13 + t^3*k1*k3*k5*k7*k8*k10*k11*k14 - t^3*k1*k4*k5*k6*k7*k8*k11*k13 - t^3*k1*k4*k5*k6*k7*k9*k11*k13 + t^3*k1*k4*k5*k7*k8*k9*k11*k14 - t^3*k2*k4*k5*k6*k7*k8*k11*k13 - t^3*k2*k4*k5*k6*k7*k9*k11*k13 + t^3*k2*k4*k5*k7*k8*k9*k11*k14 + t^3*k1*k2*k3*k6*k7*k10*k11 + t^3*k1*k2*k3*k6*k8*k10*k11 + t^3*k1*k2*k5*k6*k7*k8*k11 + t^3*k1*k2*k5*k6*k7*k9*k11 + t^3*k1*k3*k5*k6*k7*k8*k11 + t^3*k1*k3*k5*k6*k7*k9*k11 + t^3*k1*k4*k5*k6*k7*k8*k9 + t^3*k2*k4*k5*k6*k7*k8*k9 - t^2*k1*k2*k3*k6*k7*k10*k11*k13 - t^2*k1*k2*k3*k6*k8*k10*k11*k13 + t^2*k1*k2*k3*k7*k8*k10*k11*k14 - t^2*k1*k2*k5*k6*k7*k8*k11*k13 - t^2*k1*k2*k5*k6*k7*k9*k11*k13 + t^2*k1*k2*k5*k7*k8*k9*k11*k14 - t^2*k1*k3*k5*k6*k7*k8*k11*k13 - t^2*k1*k3*k5*k6*k7*k9*k11*k13 + t^2*k1*k3*k5*k7*k8*k9*k11*k14 - t^2*k1*k4*k5*k6*k7*k8*k9*k13 - t^2*k2*k4*k5*k6*k7*k8*k9*k13 + t^2*k1*k2*k3*k6*k7*k8*k11 + t^2*k1*k2*k3*k6*k7*k9*k11 + t^2*k1*k2*k5*k6*k7*k8*k9 + t^2*k1*k3*k5*k6*k7*k8*k9 - t*k1*k2*k3*k6*k7*k8*k11*k13 - t*k1*k2*k3*k6*k7*k9*k11*k13 + t*k1*k2*k3*k7*k8*k9*k11*k14 - t*k1*k2*k5*k6*k7*k8*k9*k13 - t*k1*k3*k5*k6*k7*k8*k9*k13 + t*k1*k2*k3*k6*k7*k8*k9 - k1*k2*k3*k6*k7*k8*k9*k13)/(t*k1*k2*k5*(t^3*k4*k7*k10*k11 + t^3*k4*k8*k10*k11 + t^2*k3*k7*k10*k11 + t^2*k3*k8*k10*k11 + t^2*k4*k7*k8*k11 + t^2*k4*k7*k9*k11 + t*k3*k7*k8*k11 + t*k3*k7*k9*k11 + t*k4*k7*k8*k9 + k3*k7*k8*k9))
end
@fastmath function J(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14)
    return abs(5*(-k1*k4*k5*k6*k7*k10*k11 - k1*k4*k5*k6*k8*k10*k11 - k2*k4*k5*k6*k7*k10*k11 - k2*k4*k5*k6*k8*k10*k11)*t^4 + 4*(-k1*k2*k4*k5*k7*k10*k11*k12 - k1*k2*k4*k5*k8*k10*k11*k12 + k1*k4*k5*k6*k7*k10*k11*k13 + k1*k4*k5*k6*k8*k10*k11*k13 - k1*k4*k5*k7*k8*k10*k11*k14 + k2*k4*k5*k6*k7*k10*k11*k13 + k2*k4*k5*k6*k8*k10*k11*k13 - k2*k4*k5*k7*k8*k10*k11*k14 - k1*k2*k5*k6*k7*k10*k11 - k1*k2*k5*k6*k8*k10*k11 - k1*k3*k5*k6*k7*k10*k11 - k1*k3*k5*k6*k8*k10*k11 - k1*k4*k5*k6*k7*k8*k11 - k1*k4*k5*k6*k7*k9*k11 - k2*k4*k5*k6*k7*k8*k11 - k2*k4*k5*k6*k7*k9*k11)*t^3 + 3*(-k1*k2*k3*k5*k7*k10*k11*k12 - k1*k2*k3*k5*k8*k10*k11*k12 - k1*k2*k4*k5*k7*k8*k11*k12 - k1*k2*k4*k5*k7*k9*k11*k12 + k1*k2*k5*k6*k7*k10*k11*k13 + k1*k2*k5*k6*k8*k10*k11*k13 - k1*k2*k5*k7*k8*k10*k11*k14 + k1*k3*k5*k6*k7*k10*k11*k13 + k1*k3*k5*k6*k8*k10*k11*k13 - k1*k3*k5*k7*k8*k10*k11*k14 + k1*k4*k5*k6*k7*k8*k11*k13 + k1*k4*k5*k6*k7*k9*k11*k13 - k1*k4*k5*k7*k8*k9*k11*k14 + k2*k4*k5*k6*k7*k8*k11*k13 + k2*k4*k5*k6*k7*k9*k11*k13 - k2*k4*k5*k7*k8*k9*k11*k14 - k1*k2*k3*k6*k7*k10*k11 - k1*k2*k3*k6*k8*k10*k11 - k1*k2*k5*k6*k7*k8*k11 - k1*k2*k5*k6*k7*k9*k11 - k1*k3*k5*k6*k7*k8*k11 - k1*k3*k5*k6*k7*k9*k11 - k1*k4*k5*k6*k7*k8*k9 - k2*k4*k5*k6*k7*k8*k9)*t^2 + 2*(-k1*k2*k3*k5*k7*k8*k11*k12 - k1*k2*k3*k5*k7*k9*k11*k12 + k1*k2*k3*k6*k7*k10*k11*k13 + k1*k2*k3*k6*k8*k10*k11*k13 - k1*k2*k3*k7*k8*k10*k11*k14 - k1*k2*k4*k5*k7*k8*k9*k12 + k1*k2*k5*k6*k7*k8*k11*k13 + k1*k2*k5*k6*k7*k9*k11*k13 - k1*k2*k5*k7*k8*k9*k11*k14 + k1*k3*k5*k6*k7*k8*k11*k13 + k1*k3*k5*k6*k7*k9*k11*k13 - k1*k3*k5*k7*k8*k9*k11*k14 + k1*k4*k5*k6*k7*k8*k9*k13 + k2*k4*k5*k6*k7*k8*k9*k13 - k1*k2*k3*k6*k7*k8*k11 - k1*k2*k3*k6*k7*k9*k11 - k1*k2*k5*k6*k7*k8*k9 - k1*k3*k5*k6*k7*k8*k9)*t - k1*k2*k3*k5*k7*k8*k9*k12 + k1*k2*k3*k6*k7*k8*k11*k13 + k1*k2*k3*k6*k7*k9*k11*k13 - k1*k2*k3*k7*k8*k9*k11*k14 + k1*k2*k5*k6*k7*k8*k9*k13 + k1*k3*k5*k6*k7*k8*k9*k13 - k1*k2*k3*k6*k7*k8*k9)
end

@fastmath function IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)
    return J(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,g(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k13,k14),k13,k14)*chi(kk12[1],kk12[2],g(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k13,k14))/(t*k1*k2*k5*(t^3*k4*k7*k10*k11 + t^3*k4*k8*k10*k11 + t^2*k3*k7*k10*k11 + t^2*k3*k8*k10*k11 + t^2*k4*k7*k8*k11 + t^2*k4*k7*k9*k11 + t*k3*k7*k8*k11 + t*k3*k7*k9*k11 + t*k4*k7*k8*k9 + k3*k7*k8*k9))
end

# MCKR simple. Input:Box, sample size. Output: I-hat.
@fastmath function sumo(B,NN)
    kk1=B[1]
    kk2=B[2]
    kk3=B[3]
    kk4=B[4]
    kk5=B[5]
    kk6=B[6]
    kk7=B[7]
    kk8=B[8]
    kk9=B[9]
    kk10=B[10]
    kk11=B[11]
    kk12=B[12]
    kk13=B[13]
    kk14=B[14]
    Ians=0
    COEFF=kk12[2]/(kk12[2]-kk12[1])
    t=samplor(0,kk12[2])
    k1=samplor(kk1[1],kk1[2])
    k2=samplor(kk2[1],kk2[2])
    k3=samplor(kk3[1],kk3[2])
    k4=samplor(kk4[1],kk4[2])
    k5=samplor(kk5[1],kk5[2])
    k6=samplor(kk6[1],kk6[2])
    k7=samplor(kk7[1],kk7[2])
    k8=samplor(kk8[1],kk8[2])
    k9=samplor(kk9[1],kk9[2])
    k10=samplor(kk10[1],kk10[2])
    k11=samplor(kk11[1],kk11[2])
    k13=samplor(kk13[1],kk13[2])
    k14=samplor(kk14[1],kk14[2])
    Ians+=IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)
    Ians=COEFF*Ians
    @inbounds for i=2:NN
        t=samplor(0,kk12[2])
        k1=samplor(kk1[1],kk1[2])
        k2=samplor(kk2[1],kk2[2])
        k3=samplor(kk3[1],kk3[2])
        k4=samplor(kk4[1],kk4[2])
        k5=samplor(kk5[1],kk5[2])
        k6=samplor(kk6[1],kk6[2])
        k7=samplor(kk7[1],kk7[2])
        k8=samplor(kk8[1],kk8[2])
        k9=samplor(kk9[1],kk9[2])
        k10=samplor(kk10[1],kk10[2])
        k11=samplor(kk11[1],kk11[2])
        k13=samplor(kk13[1],kk13[2])
        k14=samplor(kk14[1],kk14[2])
        delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
        Ians+=delto/i
    end
    return Ians
end

# MCKR simple. Input:Box, sample size. Output: I-hat, e-hat.
@fastmath function sumo_with_S(B,NN)
    kk1=B[1]
    kk2=B[2]
    kk3=B[3]
    kk4=B[4]
    kk5=B[5]
    kk6=B[6]
    kk7=B[7]
    kk8=B[8]
    kk9=B[9]
    kk10=B[10]
    kk11=B[11]
    kk12=B[12]
    kk13=B[13]
    kk14=B[14]
    S=0
    Ians=0
    COEFF=kk12[2]/(kk12[2]-kk12[1])
    t=samplor(0,kk12[2])
    k1=samplor(kk1[1],kk1[2])
    k2=samplor(kk2[1],kk2[2])
    k3=samplor(kk3[1],kk3[2])
    k4=samplor(kk4[1],kk4[2])
    k5=samplor(kk5[1],kk5[2])
    k6=samplor(kk6[1],kk6[2])
    k7=samplor(kk7[1],kk7[2])
    k8=samplor(kk8[1],kk8[2])
    k9=samplor(kk9[1],kk9[2])
    k10=samplor(kk10[1],kk10[2])
    k11=samplor(kk11[1],kk11[2])
    k13=samplor(kk13[1],kk13[2])
    k14=samplor(kk14[1],kk14[2])
    Ians+=IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)
    Ians=COEFF*Ians
    @inbounds for i=2:NN
        t=samplor(0,kk12[2])
        k1=samplor(kk1[1],kk1[2])
        k2=samplor(kk2[1],kk2[2])
        k3=samplor(kk3[1],kk3[2])
        k4=samplor(kk4[1],kk4[2])
        k5=samplor(kk5[1],kk5[2])
        k6=samplor(kk6[1],kk6[2])
        k7=samplor(kk7[1],kk7[2])
        k8=samplor(kk8[1],kk8[2])
        k9=samplor(kk9[1],kk9[2])
        k10=samplor(kk10[1],kk10[2])
        k11=samplor(kk11[1],kk11[2])
        k13=samplor(kk13[1],kk13[2])
        k14=samplor(kk14[1],kk14[2])
        delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
        Ians+=delto/i
        S+=delto^2*(i-1)/i
    end
    return Ians,S
end

# MCKR antithetic. Input:Box, sample size. Output: I-hat.
@fastmath function sumo_antithetic(B,NN)
    kk1=B[1]
    kk2=B[2]
    kk3=B[3]
    kk4=B[4]
    kk5=B[5]
    kk6=B[6]
    kk7=B[7]
    kk8=B[8]
    kk9=B[9]
    kk10=B[10]
    kk11=B[11]
    kk12=B[12]
    kk13=B[13]
    kk14=B[14]
    Ians=0
    centero=[(0+kk12[2])/2,(kk1[1]+kk1[2])/2,(kk2[1]+kk2[2])/2,(kk3[1]+kk3[2])/2,(kk4[1]+kk4[2])/2,(kk5[1]+kk5[2])/2,(kk6[1]+kk6[2])/2,(kk7[1]+kk7[2])/2,(kk8[1]+kk8[2])/2,(kk9[1]+kk9[2])/2,(kk10[1]+kk10[2])/2,(kk11[1]+kk11[2])/2,(kk13[1]+kk13[2])/2,(kk14[1]+kk14[2])/2]
    COEFF=kk12[2]/(kk12[2]-kk12[1])
    t=samplor(0,kk12[2])
    k1=samplor(kk1[1],kk1[2])
    k2=samplor(kk2[1],kk2[2])
    k3=samplor(kk3[1],kk3[2])
    k4=samplor(kk4[1],kk4[2])
    k5=samplor(kk5[1],kk5[2])
    k6=samplor(kk6[1],kk6[2])
    k7=samplor(kk7[1],kk7[2])
    k8=samplor(kk8[1],kk8[2])
    k9=samplor(kk9[1],kk9[2])
    k10=samplor(kk10[1],kk10[2])
    k11=samplor(kk11[1],kk11[2])
    k13=samplor(kk13[1],kk13[2])
    k14=samplor(kk14[1],kk14[2])
    Ians+=IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)
    Ians=COEFF*Ians
    t=2*centero[1]-t
    k1=2*centero[2]-k1
    k2=2*centero[3]-k2
    k3=2*centero[4]-k3
    k4=2*centero[5]-k4
    k5=2*centero[6]-k5
    k6=2*centero[7]-k6
    k7=2*centero[8]-k7
    k8=2*centero[9]-k8
    k9=2*centero[10]-k9
    k10=2*centero[11]-k10
    k11=2*centero[12]-k11
    k13=2*centero[13]-k13
    k14=2*centero[14]-k14
    delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
    Ians+=delto/2
    @inbounds for i=2:floor(NN/2)
        t=samplor(0,kk12[2])
        k1=samplor(kk1[1],kk1[2])
        k2=samplor(kk2[1],kk2[2])
        k3=samplor(kk3[1],kk3[2])
        k4=samplor(kk4[1],kk4[2])
        k5=samplor(kk5[1],kk5[2])
        k6=samplor(kk6[1],kk6[2])
        k7=samplor(kk7[1],kk7[2])
        k8=samplor(kk8[1],kk8[2])
        k9=samplor(kk9[1],kk9[2])
        k10=samplor(kk10[1],kk10[2])
        k11=samplor(kk11[1],kk11[2])
        k13=samplor(kk13[1],kk13[2])
        k14=samplor(kk14[1],kk14[2])
        delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
        Ians+=delto/(2*i+1)
        t=2*centero[1]-t
        k1=2*centero[2]-k1
        k2=2*centero[3]-k2
        k3=2*centero[4]-k3
        k4=2*centero[5]-k4
        k5=2*centero[6]-k5
        k6=2*centero[7]-k6
        k7=2*centero[8]-k7
        k8=2*centero[9]-k8
        k9=2*centero[10]-k9
        k10=2*centero[11]-k10
        k11=2*centero[12]-k11
        k13=2*centero[13]-k13
        k14=2*centero[14]-k14
        delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
        Ians+=delto/(2*i+2)
    end
    return Ians
end

# MCKR antithetic. Input:Box, sample size. Output: I-hat, e-hat.
@fastmath function sumo_antithetic_with_S(B,NN)
    kk1=B[1]
    kk2=B[2]
    kk3=B[3]
    kk4=B[4]
    kk5=B[5]
    kk6=B[6]
    kk7=B[7]
    kk8=B[8]
    kk9=B[9]
    kk10=B[10]
    kk11=B[11]
    kk12=B[12]
    kk13=B[13]
    kk14=B[14]
    S=0
    Ians=0
    centero=[(0+kk12[2])/2,(kk1[1]+kk1[2])/2,(kk2[1]+kk2[2])/2,(kk3[1]+kk3[2])/2,(kk4[1]+kk4[2])/2,(kk5[1]+kk5[2])/2,(kk6[1]+kk6[2])/2,(kk7[1]+kk7[2])/2,(kk8[1]+kk8[2])/2,(kk9[1]+kk9[2])/2,(kk10[1]+kk10[2])/2,(kk11[1]+kk11[2])/2,(kk13[1]+kk13[2])/2,(kk14[1]+kk14[2])/2]
    COEFF=kk12[2]/(kk12[2]-kk12[1])
    t=samplor(0,kk12[2])
    k1=samplor(kk1[1],kk1[2])
    k2=samplor(kk2[1],kk2[2])
    k3=samplor(kk3[1],kk3[2])
    k4=samplor(kk4[1],kk4[2])
    k5=samplor(kk5[1],kk5[2])
    k6=samplor(kk6[1],kk6[2])
    k7=samplor(kk7[1],kk7[2])
    k8=samplor(kk8[1],kk8[2])
    k9=samplor(kk9[1],kk9[2])
    k10=samplor(kk10[1],kk10[2])
    k11=samplor(kk11[1],kk11[2])
    k13=samplor(kk13[1],kk13[2])
    k14=samplor(kk14[1],kk14[2])
    Ians+=IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)
    Ians=COEFF*Ians
    t=2*centero[1]-t
    k1=2*centero[2]-k1
    k2=2*centero[3]-k2
    k3=2*centero[4]-k3
    k4=2*centero[5]-k4
    k5=2*centero[6]-k5
    k6=2*centero[7]-k6
    k7=2*centero[8]-k7
    k8=2*centero[9]-k8
    k9=2*centero[10]-k9
    k10=2*centero[11]-k10
    k11=2*centero[12]-k11
    k13=2*centero[13]-k13
    k14=2*centero[14]-k14
    delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
    Ians+=delto/2
    S+=delto^2/2
    @inbounds for i=2:floor(NN/2)
        t=samplor(0,kk12[2])
        k1=samplor(kk1[1],kk1[2])
        k2=samplor(kk2[1],kk2[2])
        k3=samplor(kk3[1],kk3[2])
        k4=samplor(kk4[1],kk4[2])
        k5=samplor(kk5[1],kk5[2])
        k6=samplor(kk6[1],kk6[2])
        k7=samplor(kk7[1],kk7[2])
        k8=samplor(kk8[1],kk8[2])
        k9=samplor(kk9[1],kk9[2])
        k10=samplor(kk10[1],kk10[2])
        k11=samplor(kk11[1],kk11[2])
        k13=samplor(kk13[1],kk13[2])
        k14=samplor(kk14[1],kk14[2])
        delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
        Ians+=delto/(2*i+1)
        S+=(delto^2)*2*i/(2*i+1)
        t=2*centero[1]-t
        k1=2*centero[2]-k1
        k2=2*centero[3]-k2
        k3=2*centero[4]-k3
        k4=2*centero[5]-k4
        k5=2*centero[6]-k5
        k6=2*centero[7]-k6
        k7=2*centero[8]-k7
        k8=2*centero[9]-k8
        k9=2*centero[10]-k9
        k10=2*centero[11]-k10
        k11=2*centero[12]-k11
        k13=2*centero[13]-k13
        k14=2*centero[14]-k14
        delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
        Ians+=delto/(2*i+2)
        S+=(delto^2)*(2*i+1)/(2*i+2)
    end
    return Ians,S
end
