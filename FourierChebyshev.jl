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

function plot_fouriercoeffs(bc)
    N = (length(bc)+1)/2 # 2N-1
    plot(-N+1:N-1,abs.(bc),yscale=:log10,
    legend=false,
    xlabel = "\$k\$",
    ylabel = "\$|\\bar{c}_k\\,|\$",
    title = "Absolute values of Fourier coefficients"
    )
end