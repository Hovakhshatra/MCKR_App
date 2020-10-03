# The Simmple and Antithetic MCKR for example number 3 with parameters equipped with Uniform distibution.

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

@fastmath function g2(t, k1)
    return 25 / 4 - 100 * t^5 + 100 * t^4 * k1 + 450 * t^4 - 450 * t^3 * k1 - 525 * t^3 + 575 * t^2 * k1 - 75 / 2 * t^2 - 375 / 2 * k1 * t + 575 / 2 * t;
end
@fastmath function J(t, k1)
    return abs(5 * t^4 + 4 * (-k1 - 9 / 2) * t^3 + 3 * ((9 * k1) / 2 + 21 / 4) * t^2 + 2 * (-(23 * k1) / 4 + 3 / 8) * t + (15 * k1) / 8 - 23 / 8)
end

@fastmath function SummandWithoutCOEFF(t, k1, ak2, bk2)
    return 100 * J(t, k1) * chi(ak2, bk2, g2(t, k1))
end

# MCKR simple. Input:Box, sample size. Output: I-hat.
@fastmath function sumo(B, NN)
    kk1=B[1]
    kk2=B[2]
    Ians = 0
    COEFF = 1 / (kk2[2] - kk2[1])
    t = samplor(0, 1)
    k1 = samplor(kk1[1], kk1[2])
    Ians += SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2
    Ians = COEFF * Ians
    @inbounds for i = 2:NN
        t = samplor(0, 1)
        k1 = samplor(kk1[1], kk1[2])
        delto = COEFF * ( SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) ) / t^2 - Ians
        Ians += delto / i
    end
    return Ians
end

# MCKR simple. Input:Box, sample size. Output: I-hat, e-hat.
@fastmath function sumo_with_S(B, NN)
    kk1=B[1]
    kk2=B[2]
    S = 0
    Ians = 0
    COEFF = 1 / (kk2[2] - kk2[1])
    t = samplor(0, 1)
    k1 = samplor(kk1[1], kk1[2])
    Ians += SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2
    Ians = COEFF * Ians
    @inbounds for i = 2:NN
        t = samplor(0, 1)
        k1 = samplor(kk1[1], kk1[2])
        delto = COEFF * ( SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) ) / t^2 - Ians
        Ians += delto / i
        S += delto^2 * (i - 1) / i
    end
    return Ians, S
end

# MCKR antithetic. Input:Box, sample size. Output: I-hat.
@fastmath function sumo_antithetic_with_S(B, NN)
    kk1=B[1]
    kk2=B[2]
    Ians = 0
    centero = [1 / 2,(kk1[1] + kk1[2]) / 2]
    COEFF = 1 / (kk2[2] - kk2[1])
    t = samplor(0, 1)
    k1 = samplor(kk1[1], kk1[2])
    Ians += SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2
    Ians = COEFF * Ians
    t = 2 * centero[1] - t
    k1 = 2 * centero[2] - k1
    delto = SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2 - Ians
    Ians += delto / 2
    @inbounds for i = 1:floor(NN / 2)
        t = samplor(0, 1)
        k1 = samplor(kk1[1], kk1[2])
        delto = SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2 - Ians
        Ians += delto / (2 * i + 1)
        t = 2 * centero[1] - t
        k1 = 2 * centero[2] - k1
        delto = SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2 - Ians
        Ians += delto / (2 * i + 2)
    end
    return Ians
end

# MCKR antithetic. Input:Box, sample size. Output: I-hat, e-hat.
@fastmath function sumo_antithetic_with_S(B, NN)
    kk1=B[1]
    kk2=B[2]
    S = 0
    Ians = 0
    centero = [1 / 2,(kk1[1] + kk1[2]) / 2]
    COEFF = 1 / (kk2[2] - kk2[1])
    t = samplor(0, 1)
    k1 = samplor(kk1[1], kk1[2])
    Ians += SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2
    Ians = COEFF * Ians
    t = 2 * centero[1] - t
    k1 = 2 * centero[2] - k1
    delto = SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2 - Ians
    Ians += delto / 2
    S += delto^2 / 2
    @inbounds for i = 1:floor(NN / 2)
        t = samplor(0, 1)
        k1 = samplor(kk1[1], kk1[2])
        delto = SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2 - Ians
        Ians += delto / (2 * i + 1)
        S += 2 * i * delto^2 / (2 * i + 1)
        t = 2 * centero[1] - t
        k1 = 2 * centero[2] - k1
        delto = SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2 - Ians
        Ians += delto / (2 * i + 2)
        S += (2 * i + 1) * delto^2 / (2 * i + 2)
    end
    return Ians, S
end
