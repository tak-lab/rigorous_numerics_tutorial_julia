### Fourier functions
using FFTW, Plots
function fouriercoeffs(f, N, I=[0,2π])
    # f: any periodic function on I
    # N: size of Fourier series
    a = I[1]; b = I[2]
    h = (b-a)/(2N-1)
    j = 0:2N-2
    xⱼ = a .+ j*h
    fⱼ = f.(xⱼ);
    return fftshift(fft(fⱼ))/(2N-1)
end

function odefouriercoeffs(f, N, I, n=1)
    a = I[1]; b = I[2];
    h = (b-a)/(2N-1)
    j = 0:2N-2
    xⱼ = a .+ j*h
    fⱼ = f(xⱼ)[n,:]
    return fftshift(fft(fⱼ))/(2N-1)
end

function plot_fourier(bc, I=[0,2π])
    # bc: Fourier coefficients
    a = I[1]; b = I[2]
    N = (length(bc)+1)/2 # 2N-1
    n_pad = 500
    bc_pad = [zeros(n_pad);bc;zeros(n_pad)]
    N_pad = N + n_pad
    h_pad = (b-a)/(2N_pad-1)
    xj_pad = a .+ h_pad*(0:2N_pad-2)
　　fNj_pad = real((2N_pad-1)*ifft(ifftshift(bc_pad)))
　　plot(xj_pad, fNj_pad, legend=false, xlabel = "\$x\$", ylabel = "\$f(x)\$")
end

function plot_fourier!(bc, I=[0,2π];label="")
    # bc: Fourier coefficients
    a = I[1]; b = I[2]
    N = (length(bc)+1)/2 # 2N-1
    n_pad = 500
    bc_pad = [zeros(n_pad);bc;zeros(n_pad)]
    N_pad = N + n_pad
    h_pad = (b-a)/(2N_pad-1)
    xj_pad = a .+ h_pad*(0:2N_pad-2)
　　fNj_pad = real((2N_pad-1)*ifft(ifftshift(bc_pad)))
　　plot!(xj_pad, fNj_pad, label=label)
end

function plot_fouriercoeffs(bc)
    N = (length(bc)+1)/2 # 2N-1
    plot(-N+1:N-1,abs.(bc),yscale=:log10,
    legend=false,
    xlabel = "\$k\$",
    ylabel = "\$|\\bar{c}_k\\,|\$",
    title = "Absolute values of Fourier coefficients"
    )
end

function plot_solution(u, index) # u = [ω, a_{-N+1}, ..., a_0, ..., a_{N-1}], length(u) = 2N
    # index = 1: profile of solution
    #         2: Fourier mode
    #         3: phase profile
    ω = real(u[1])
    L = 2π / ω
    a = u[2:end]
    N = length(u)/2 # N: size of Fourier
    n_pad = 500
    a_pad = [zeros(n_pad);a;zeros(n_pad)]
    N_pad = N + n_pad    
    dx = L/(2*N_pad-1)
    x = dx*(0:2*N_pad-2)
    if index == 1
    # Plot profile:
        plot(plot_fourier(a, [0,L]),
        # plot(x,real((2N_pad-1)*ifft(ifftshift(a_pad))),
            xlabel = "\$t\$",
            ylabel = "\$x\\,(t)\$",
            line   = 1.6,
            title  = "Profile of solution",
            size   = (720,400),
            legend = false,
        )
    elseif index == 2
    # Plot Fourier coefficients:
        plot(plot_fouriercoeffs(a),
        # plot((-N+1):(N-1),abs.(a),yscale=:log10,
            xlabel = "\$k\$",
            ylabel = "\$|a_k\\,|\$",
            line   = 1.6,
            title  = "Absolute values of Fourier coefficients",
            size   = (720,400),
            legend = false,
        )
    elseif index == 3
    # Plot phase:
      k = (-N_pad+1):(N_pad-1)
      plot(real((2N_pad-1)*ifft(ifftshift(a_pad))),real((2N_pad-1)*ifft(ifftshift(a_pad.*(im*k*ω)))),
            xlabel = "\$x(t)\$",
            ylabel = "\$\\dot{x}\\,(t)\$",
            line   = 1.6,
            title  = "Phase plot of a numerical solution",
            size   = (720,400),
            legend = false,
        )
    end
end

function plot_solution!(u)
    L = 2π/real(u[1])
    a = u[2:end]
    N = length(u)/2
    n_pad = 1000
    a_pad = [zeros(n_pad);a;zeros(n_pad)]
    N_pad = N+n_pad
    k = (-N_pad+1):(N_pad-1)
    dx = L/(2*N_pad-1)
    x = dx*(0:2*N_pad-2)
    plot!(real((2N_pad-1)*ifft(ifftshift(a_pad))),real((2N_pad-1)*ifft(ifftshift(a_pad.*(im*k)))),line=1.6,)
end

function powerconvfourier(a::Vector{Complex{T}},p) where T
    M = Int((length(a)+1)/2)
    N = (p-1)*M
    ta = [zeros(N,1);a;zeros(N,1)] # 1. Padding zeros: size(ta) = 2pM-1
    tb = ifft(ifftshift(ta)) # 2. IFFT of ta
    tbᵖ = tb.^p # 3. tb*^tb
    cᵖ = fftshift(fft(tbᵖ))*(2.0*p*M-1)^(p-1)
    return cᵖ[N+1:end-N], cᵖ[p:end-(p-1)]# return (truncated, full) version
end


### Chebyshev functions
function chebpts(n, a=-1, b=1) # n: maximum order of Chebyshev polynomials
    tt = range(0, stop=π, length=n+1)
    x = cos.(tt)
    return (1.0 .- x).*a/2 + (1.0 .+ x).*b/2
end

function chebcoeffs(f,M,I=[-1,1])
    a = I[1]; b = I[2]
    n = M-1
    cpts  = chebpts(n, a, b)
    fvals = f.(cpts)
    FourierCoeffs = real(fft([fvals;reverse(fvals[2:end-1])]))
    ChebCoeffs = FourierCoeffs[1:n+1]/n
    ChebCoeffs[1] = ChebCoeffs[1]/2
    ChebCoeffs[end] = ChebCoeffs[end]/2
    return ChebCoeffs # return Two-sided Chebyshev
end

function cheb(f,I=[-1;1]; tol = 5e-15,Nmax = 10000)
    a = I[1]; b = I[2]; m = 0.5*(a+b); r = 0.5*(b-a); x = rand(5)
    x1 = m .+ x*r; x2 = m .- x*r
    if f.(x1) ≈ f.(x2)
        odd_even = 1 # even function: 1
    elseif f.(x1) ≈ -f.(x2)
        odd_even = -1 #  odd function: -1
    else
        odd_even = 0 # otherwise: 0
    end
    i = 3
    schbc = 0 # sampling chebyshev coefficients
    while true
        schbc = chebcoeffs(f,2^i+1,I)
        if all(abs.(schbc[end-2:end]) .< tol) || (2^i+1 > Nmax) 
            break
        end
        i += 1
    end    
    M = findlast(abs.(schbc) .> tol)
    cc = chebcoeffs(f,M,I)
    if odd_even == 1 # even function
        cc[2:2:end] .= 0
    elseif odd_even == -1 # odd function
        cc[1:2:end] .= 0
    end
    return cc # return Two-sided Chebyshev
end

function plot_chebcoeffs(f)
    zero_ind = findall(x->x==0, f)
    f[zero_ind] .= f[zero_ind .+ 1]
    plot(0:length(f)-1, abs.(f),
        yscale=:log10,
        title="Chebyshev coefficients",
        xlabel="Degree of Chebyshev polynomial",
        ylabel="Magnitude of coefficient",
        size       = (800,400),
        legend     = false,
    )
end

function clenshaw(a,x) # Clenshaw's algorithm
# a: (Two-sided) Chbyshev coefficients
# x: evaluating points in [-1,1]
    n = length(a)-1
    bk1 = 0.0
    bk2 = 0.0
    x = 2x
    for r = (n+1):-2:3
        bk2 = x.*bk1 .- bk2 .+ a[r]
        bk1 = x.*bk2 .- bk1 .+ a[r-1]
    end
    if isodd(n)        
        b2 = x.*bk1 .- bk2 .+ a[2]
        bk2 = bk1 # b3
        bk1 = b2
    end
    return -bk2 .+ 0.5x .* bk1 .+ a[1] # y = c(1) + .5*x.*bk1 - bk2;
end

function clenshaw_secondkind(a,x) # Clenshaw's algorithm
# a: (Two-sided) Chbyshev coefficients
# x: evaluating points in [-1,1]
    n = length(a)-1
    bk0 = 0
    bk1 = 0
    for r = (n+1):-1:1
        tmp = 2x.*bk0 .- bk1 .+ a[r]
        bk1 = bk0
        bk0 = tmp
    end
    return bk0 #.+ bk1*(-x)
end

function eval_cheb(a,x; I=[-1,1])
# a: (Two-sided) Chbyshev coefficients
# x: evaluating points in domain
    ξ = 2*(x.-I[1])/(I[2]-I[1]) .- 1
    return clenshaw(a,ξ)
end

function eval_cheb_naive(ChebCoeffs_twosided,x; I=[-1,1])
    M = length(ChebCoeffs_twosided) # M: size of chebyshev
    a = I[1]; b = I[2]
    k = 0:M-1
    ξ = 2*(x.-a)/(b-a) .- 1
    return cos.(Vector(k)' .* acos.(ξ)) * ChebCoeffs_twosided
end

function eval_cheb_bc(ChebCoeffs_twosided,x,n=200; I=[-1,1]) # Barycentric interportion formula
    M = length(ChebCoeffs_twosided) # M: size of chebyshev
    a = I[1]; b = I[2]
    k = 0:M-1
    ξⱼ = chebpts(n)
    xc = (1.0 .- ξⱼ)*a/2 + (1.0 .+ ξⱼ)*b/2 # Chebyshev points in [a,b]
    fxc = cos.(Vector(k)' .* acos.(ξⱼ)) * ChebCoeffs_twosided
    valnum = length(x)
    ξ = 2*(x.-a)/(b-a) .- 1
    # ξ = range(-1,stop=1,length=valnum)
    # x = (1.0 .- ξ)*a/2 + (1.0 .+ ξ)*b/2
    λ = [1/2; ones(n-1); 1/2] .* (-1).^(0:n)

    numer = zeros(valnum)
    denom = zeros(valnum)
    exact = zeros(Bool,valnum)

    for j = 1:n+1
        xdiff = x .- xc[j]
        temp = λ[j] ./ xdiff
        numer += temp * fxc[j]
        denom += temp
        exact[xdiff.==0] .= true
    end

    fx = numer ./ denom
    jj = findall(exact)
    fx[jj] = f.(x[jj])
    fx[jj] = cos.(Vector(k)' .* acos.(ξ[jj])) * ChebCoeffs_twosided
    return fx
end

function plot_cheb(ChebCoeffs_twosided; I=[-1,1],title="",label="",legend=true) # Input: Two-sided Chebyshev
    # M = length(ChebCoeffs_twosided) # M: size of chebyshev
    # a = I[1]; b = I[2]; 
    x = range(I[1],stop=I[2],length=5000)
    # ξ = 2*(x.-a)/(b-a) .- 1
    fx = eval_cheb(ChebCoeffs_twosided,x,I=I)
    plot(x, fx, legend=legend, label=label, title=title, xlabel="\$x\$",ylabel="\$f(x)\$")
end

function plot_cheb!(ChebCoeffs_twosided; I=[-1,1],title="",label="",legend=true) # Input: Two-sided Chebyshev
    # M = length(ChebCoeffs_twosided) # M: size of chebyshev
    # a = I[1]; b = I[2]; 
    x = range(I[1],stop=I[2],length=5000)
    # ξ = 2*(x.-a)/(b-a) .- 1
    fx = eval_cheb(ChebCoeffs_twosided,x,I=I)
    plot!(x, fx, legend=legend, label=label, title=title, xlabel="\$x\$",ylabel="\$f(x)\$")
end

function chebdiff(a; I=[-1,1])# Input is Two-sided
    M = length(a)
    b = zeros(M+1)
    for r = M-1:-1:1
        b[r] = b[r+2] + 2*r*a[r+1]
    end
    b[1] /= 2.0
    return b[setdiff(1:end,end)]*(2/(I[2]-I[1])) # Output is Two-sided
end

function chebdiff_oneside(a; I=[-1,1])# Input is One-sided
    M = length(a)
    b = zeros(M+1)
    for r = M-1:-1:1
        b[r] = b[r+2] + 2*r*a[r+1]
    end
    return b[setdiff(1:end,end)]*(2/(I[2]-I[1])) # Output is One-sided
end

function chebdiff_secondkind(a; I=[-1,1])　# Input is Two-sided
    M = length(a)
    b = zeros(M-1)
    for n = 0:M-2
        b[n+1] = (n+1)*a[n+2]
    end
    return b*(2/(I[2]-I[1]))　# Output is second kind (Two-sided)
end

function chebindefint(a; I=[-1,1])# Input is Two-sided
    M = length(a)
    a_ext = zeros(M+2)
    a_ext[1] = 2*a[1]
    a_ext[2:M] = a[2:M]
    A = zeros(M+1)
    for n = 1:M
        A[n+1] = (a_ext[n] - a_ext[n+2])/(2n)
    end
    # A[1] = sum(A[2:2:end]) - sum(A[3:2:end]) # takes the value 0 at the left endpoint
    return A * (I[2]-I[1])/2
end

function chebint(a; I=[-1,1])# Input is Two-sided
    M = length(a)
    n = 0:2:M-1
    return sum(2a[1:2:end]./(1.0 .- n.^2))*((I[2]-I[1])/2)
end