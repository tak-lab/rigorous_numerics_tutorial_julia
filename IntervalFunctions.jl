### Interval Matrix Multiplication
function ufp(P)
    u = 2.0^(-53);
    ϕ = 2.0^52 + 1;
    q = ϕ * P;
    T = (1 - u)*q;
    return abs(q - T);
end

function succ(c)
    s_min = 2.0^-1074;
    u = 2.0^-53;
    ϕ = u * (1.0 + 2.0 * u);
    if abs(c) >= (1. / 2.) * u^(-2) * s_min # 2^(-969)(Float64)
        e = ϕ * abs(c);
        succ = c + e;
    elseif abs(c) < (1. / 2.) * u^(-1) * s_min # 2^(-1022)(Float64)
        succ = c + s_min;
    else
        C = u^(-1) * c;
        e = ϕ * abs(C);
        succ = (C + e) * u;
    end
    return succ
end

function pred(c)
    s_min = 2.0^-1074;
    u = 2.0^-53;
    ϕ = u * (1.0 + 2.0 * u);
    if abs(c) >= (1. / 2.) * u^(-2) * s_min # 2^(-969)(Float64)
        e = ϕ * abs(c);
        pred = c - e;
    elseif abs(c) < (1. / 2.) * u^(-1) * s_min # 2^(-1022)(Float64)
        pred = c - s_min;
    else
        C = u^(-1) * c;
        e = ϕ * abs(C);
        pred = (C - e) * u;
    end
    return pred
end

function mm_ufp(A_mid, B_mid) # A_mid, B_mid: Point matrix
    u = 2.0^(-53);
    realmin = 2.0^(-1022);
    n = size(A_mid,2);

    if(2*(n+2)*u>=1)
        error("mm_ufp is failed!(2(n+2)u>=1)")
    end
    # C_mid = A_mid * B_mid;
    # C_rad = (n+2) * u * ufp.(abs.(A_mid)*abs.(B_mid)) .+ realmin;
    # return C_mid, C_rad;
    return A_mid * B_mid, (n+2) * u * ufp.(abs.(A_mid)*abs.(B_mid)) .+ realmin
end

function imm_ufp(A_mid, A_rad, B_mid, B_rad) # A = <A_mid, A_rad>, B = <B_mid, B_rad>: Interval matrix
    u = 2.0^(-53);
    realmin = 2.0^(-1022);
    n = size(A_mid,2);

    if(2*(n+2)*u>=1)
        error("mm_ufp is failed!(2(n+2)u>=1)")
    end
#     C, R = mm_ufp(A_mid,B_mid);
    # C_mid = A_mid * B_mid;
    R = (n+2) * u * ufp.(abs.(A_mid)*abs.(B_mid)) .+ realmin;

#     T_1, T_2 = mm_ufp(abs.(A_mid), B_rad);
    T1 = abs.(A_mid) * B_rad;
    T2 = (n+2)*u*ufp.(T1) .+ realmin;

#     T_3 = succ.(abs.(B_mid)+B_rad);
    T3 = succ.(abs.(B_mid)+B_rad)

#     T_4, T_5 = mm_ufp(A_r, T_3);
    T4 = A_rad * T3;
    T5 = (n+2)*u*ufp.(T4) .+ realmin;

    rad_sum = R + T1 + T2 + T4 + T5;

    # C_rad = succ.(rad_sum + 4*u*ufp.(rad_sum));

    # return C_mid, C_rad;
    return A_mid * B_mid, succ.(rad_sum + 4*u*ufp.(rad_sum))
end

using IntervalArithmetic
function int_mul(A::Matrix{T}, B::Matrix{T}) where T
    Cmid, Crad = mm_ufp(A, B);
    return Cmid .± Crad
end

function int_mul(A::Matrix{Interval{T}}, B::Matrix{T}) where T
    Cmid, Crad = imm_ufp(mid.(A), radius.(A), B, zeros(size(B)));
    return Cmid .± Crad
end

function int_mul(A::Matrix{T}, B::Matrix{Interval{T}}) where T
    Cmid, Crad = imm_ufp(A, zeros(size(A)), mid.(B), radius.(B));
    return Cmid .± Crad
end

function int_mul(A::Matrix{Complex{T}}, B::Matrix{T}) where T
    Ar = real.(A); Ai = imag.(A); # (Ar + im*Ai)*B = Ar*B + im*(Ai*B)
    return int_mul(Ar, B) + im * int_mul(Ar, B)
end

function int_mul(A::Matrix{T}, B::Matrix{Complex{T}}) where T
    Br = real.(B); Bi = imag.(B); # A*(Br + im*Bi) = A*Br + im*(A*Bi)
    return int_mul(A, Br) + im * int_mul(A, Bi)
end

function int_mul(A::Matrix{Complex{T}}, B::Matrix{Complex{T}}) where T
    Ar = real.(A); Ai = imag.(A); Br = real.(B); Bi = imag.(B);
    # (Ar + im*Ai)*(Br + im*Bi) = (Ar*Br - Ai*Bi) + im*(Ar*Bi + Ai*Br)
    return (int_mul(Ar,Br) - int_mul(Ai, Bi)) + im * (int_mul(Ar, Bi) + int_mul(Ai, Br))
end


### Interval Linear system solver
