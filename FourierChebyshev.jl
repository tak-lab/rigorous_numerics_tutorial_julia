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

function cheb(f,I=[-1;1];tol = 5e-15,Nmax = 10000)
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

function eval_cheb(ChebCoeffs_twosided,x,n=200)
    M = length(ChebCoeffs_twosided) # M: size of chebyshev
    a = x[1]; b = x[end]
    k = 0:M-1
    ξⱼ = chebpts(n)
    xc = (1.0 .- ξⱼ)*a/2 + (1.0 .+ ξⱼ)*b/2 # Chebyshev points in [a,b]
    fxc = cos.(Vector(k)' .* acos.(ξⱼ)) * ChebCoeffs_twosided
    valnum = length(x)
    ξ = 2*(x.-a)/(b-a) .- 1;
    # ξ = range(-1,stop=1,length=valnum)
    x = (1.0 .- ξ)*a/2 + (1.0 .+ ξ)*b/2
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

function plot_cheb(ChebCoeffs_twosided;n=200,I=[-1,1],title="",label="",legend=true) # Input: Two-sided Chebyshev
    # M = length(ChebCoeffs_twosided) # M: size of chebyshev
    a = I[1]; b = I[2]; 
    x = range(a,stop=b,length=5000)
    fx = eval_cheb(ChebCoeffs_twosided,x,n)
#     k = 0:M-1
#     ξⱼ = chebpts(n)
#     xc = (1.0 .- ξⱼ)*a/2 + (1.0 .+ ξⱼ)*b/2 # Chebyshev points in [a,b]
#     fxc = cos.(Vector(k)' .* acos.(ξⱼ)) * ChebCoeffs_twosided
    
#     ξ = range(-1,stop=1,length=2000)
#     x = (1.0 .- ξ)*a/2 + (1.0 .+ ξ)*b/2
#     λ = [1/2; ones(n-1); 1/2] .* (-1).^(0:n)

#     numer = zeros(size(x))
#     denom = zeros(size(x))
#     exact = zeros(Bool,size(x))

#     for j = 1:n+1
#         xdiff = x .- xc[j]
#         temp = λ[j] ./ xdiff
#         numer += temp * fxc[j]
#         denom += temp
#         exact[xdiff.==0] .= true
#     end

#     fx = numer ./ denom
#     jj = findall(exact)
#     fx[jj] = f.(x[jj])
    plot(x, fx, legend=legend, label=label, title=title, xlabel="\$x\$",ylabel="\$f(x)\$")
end

function plot_cheb!(ChebCoeffs_twosided;n=200,I=[-1,1],title="",label="",legend=true) # Input: Two-sided Chebyshev
    # M = length(ChebCoeffs_twosided) # M: size of chebyshev
    a = I[1]; b = I[2]; 
    x = range(a,stop=b,length=5000)
    fx = eval_cheb(ChebCoeffs_twosided,x,n)
#     k = 0:M-1
#     ξⱼ = chebpts(n)
#     xc = (1.0 .- ξⱼ)*a/2 + (1.0 .+ ξⱼ)*b/2 # Chebyshev points in [a,b]
#     fxc = cos.(Vector(k)' .* acos.(ξⱼ)) * ChebCoeffs_twosided
    
#     ξ = range(-1,stop=1,length=2000)
#     x = (1.0 .- ξ)*a/2 + (1.0 .+ ξ)*b/2
#     λ = [1/2; ones(n-1); 1/2] .* (-1).^(0:n)

#     numer = zeros(size(x))
#     denom = zeros(size(x))
#     exact = zeros(Bool,size(x))

#     for j = 1:n+1
#         xdiff = x .- xc[j]
#         temp = λ[j] ./ xdiff
#         numer += temp * fxc[j]
#         denom += temp
#         exact[xdiff.==0] .= true
#     end

#     fx = numer ./ denom
#     jj = findall(exact)
#     fx[jj] = f.(x[jj])
    plot!(x, fx, legend=legend, label=label, title=title, xlabel="\$x\$",ylabel="\$f(x)\$")
end