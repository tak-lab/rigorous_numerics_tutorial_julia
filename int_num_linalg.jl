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

function mm_ufp(A_mid, B_mid)
    u = 2.0^(-53);
    realmin = 2.0^(-1022);
    n = size(A_mid,2);

    if(2*(n+2)*u>=1)
        error("mm_ufp is failed!(2(n+2)u>=1)")
    end
    C_mid = A_mid * B_mid;
    C_rad = (n+2) * u * ufp.(abs.(A_mid)*abs.(B_mid)) .+ realmin;

    return C_mid, C_rad;
end

function imm_ufp(A_mid, A_rad, B_mid, B_rad)
    u = 2.0^(-53);
    realmin = 2.0^(-1022);
    n = size(A_mid,2);

    if(2*(n+2)*u>=1)
        error("mm_ufp is failed!(2(n+2)u>=1)")
    end
#     C, R = mm_ufp(A_mid,B_mid);
    C_mid = A_mid * B_mid;
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

    C_rad = succ.(rad_sum + 4*u*ufp.(rad_sum));

    return C_mid, C_rad;
end
