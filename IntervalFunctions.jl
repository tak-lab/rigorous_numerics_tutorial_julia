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

# USE IntervalArithmetic.jl
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

function int_mul(A::Matrix{Interval{T}}, B::Matrix{Interval{T}}) where T
    Cmid, Crad = imm_ufp(mid.(A), radius.(A), mid.(B), radius.(B));
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



### Verify FFT using Interval Arithmetic
function verifyfft(z::Vector{T}, sign=1) where T
    n = length(z); col = 1; array1 = true
    if n==1
        Z = map(T,z)
        return Z
    else
        isrow_ = false
    end
    log2n = Int(round(log2(n))) #check dimension
    if 2^log2n ≠ n # return error if n is not the powers of 2
        error("length must be power of 2")
    end
    #bit-reversal
    f = 2^(log2n-1)
    v = [0;f]
    for k = 1:log2n-1
        f = f >> 1
        v = append!(v,f.+v)
    end
    z2 = zeros(n,col)
    if isa(real(z[1]),Interval)
        z2 = map(T,z2)
    end
    # replace z
    for j = 1: n
        z2[j,:] = z[v[j]+1,:]
    end
    #Danielson-Lanczos algorithm
    Z = complex(map(Interval,z2))
    Index = reshape([1:n*col;],n,col)

    theta = map(Interval,sign * (0:n-1)/n); # division exact because n is power of 2
    Phi = cospi.(theta) + im*sinpi.(theta) # SLOW?

    v = [1:2:n;]
    w = [2:2:n;]
    t = Z[w,:]
    Z[w,:]  = Z[v,:] - t
    Z[v,:]  = Z[v,:] + t
    for index　in 1: (log2n-1)    
        m = 2^index
        m2 = 2*m
        vw = reshape([1:n;],m2,Int(n/m2))
        v = vw[1: m, :]
        w = vw[m+1: m2, : ]
        indexv = reshape(Index[v[:],:],m,Int(col*n/m2))
        indexw = reshape(Index[w[:],:],m,Int(col*n/m2))
        Phi1 = repeat(Phi[1:Int(n/m):end],outer=[1,Int(col*n/m2)])
        t = Phi1 .*  Z[indexw]
        Z[indexw] = Z[indexv] - t 
        Z[indexv] = Z[indexv] + t
    end
    reverse(Z[2:end,:],dims=2)
     if sign==-1
        Z = Z/n
    end
    if isrow_
        Z = transpose(Z)　#transpose of Z
    end
    if array1
        Z = Z[:,1]
    end
    return Z
end

### Rigorous convolution algorithm via FFT
function powerconvfourier(a::Vector{Complex{Interval{T}}},p) where T
    M = Int((length(a)+1)/2) # length(a) = 2M-1
    N = (p-1)*M
    ia = map(Interval, a)

    length_ia = 2*p*M-1
    length_ia_ext = nextpow(2,length_ia)# 2pM-2+2L
    
    L = Int((length_ia_ext - length_ia + 1)/2)
    
    # step.1 : padding (p-1)M + L zeros for each sides
    ia_ext = map(Complex{Interval},zeros(length_ia_ext))
    ia_ext[L+N+1:end-L-N+1] = ia  #\tilda{a}

    # step.2 : inverse fft
    ib_ext = verifyfft(ifftshift(ia_ext), -1) #sign = -1 : ifft
    
    # step.3 : power p elementwisely
    ib_extᵖ = ib_ext.^p
    
    # step.4 : fft with rescaling
    ic_extᵖ = fftshift(verifyfft(ib_extᵖ, 1)) * length_ia_ext^(p-1)  #sign = 1 : fft
    
#     return ic_extᵖ,ic_extᵖ
    return ic_extᵖ[L+N+1:end-N-L+1], ic_extᵖ[L+p:end-(L+p-2)] # return (truncated, full) version
end