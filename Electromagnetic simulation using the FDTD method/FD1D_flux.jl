#=
    1D FDTD simulation
    New Formulation using flux density
    The Fourier Transform has been added
    In frequency dependent material
=#

using PyPlot
const c = 299792458.0  # speed of light [m/s]
const μ0 = 4e-7pi      # vacuum permeability
const ε0 = 1/(μ0*c^2) # vacuum permittivity
const εr = 2.0        # relative permittivity of dielectric
const σ = 0.01         # conductivity of dielectric
const KE = 200         # KE is the number of cells to be used
const τ = 0.001e-6
const χ = 2.0
const kc = 5           # input pulse location
const Nfreq = 5        # number of frequency in the interest

## variables for solution
ga = ones(KE)
gb = zeros(KE)
gc = zeros(KE)
ex = zeros(KE)
dx = zeros(KE)
hy = zeros(KE)
ix = zeros(KE)
sx = zeros(KE)

## parameters for gaussian pulse
t0 = 50.0     # Center of the incident pulse (gaussian pulse)
spread = 10.0 # Width of the incident pulse  

## parameters for sinusoidal pulse
freq_in = 700e6       # frequency [Hz]
λ = c/freq_in      # wavelength [m]
Nsampling = 10      # number of points per wavelength

ddx = 0.01
#ddx = λ/Nsampling  # cell size
dt = ddx/2c         # time step

del_exp = exp(-dt/τ)

## variables for DFT
real_pt = zeros(Nfreq, KE)   # Real and imaginary parts
imag_pt = zeros(Nfreq, KE)   # of the Fourier transform
ampn = zeros(Nfreq, KE)      # Amplitude and phase of the
phasen = zeros(Nfreq, KE)    # Fourier transform
real_in = zeros(Nfreq)       # Fourier trans. of input pulse
imag_in = zeros(Nfreq)
amp_in = zeros(Nfreq)
phase_in = zeros(Nfreq)
freq = zeros(Nfreq)
freq[1] = 50e6
freq[2] = 200e6
freq[3] = 500e6
arg = [2pi*freq[i]*dt for i in 1:Nfreq]

T = 0
NSTEPS = 1

## parameters for absorbing boundary condition
ex_low_m1 = ex_low_m2 = ex_high_m1 = ex_high_m2 = 0.0

## filling parameters in dielectric area
dielectric_range = 0.5KE:KE # range of dielctric
ga[dielectric_range] = 1.0/(εr + σ*dt/ε0 + χ*dt/τ)
gb[dielectric_range] = σ*dt/ε0
gc[dielectric_range] = χ*dt/τ


## Main part of the program
while NSTEPS > 0
    @printf("nsteps --> ")
    NSTEPS = int(readline())
    @printf("T = %d \n", T+NSTEPS)
    NSTEPS >= 1 || (break)

    for n in 1:NSTEPS
        T += 1

        ## Calculate the Dx field
        for k in 2:KE
            dx[k] += 0.5(hy[k-1] - hy[k])
        end

        ## Put a Gaussian pulse at the low end
        pulse = exp(-0.5((t0-T)/spread)^2.0)
        ## Sinusoidal pulse
        # pulse = sin(2pi*freq_in*dt*T)
        dx[kc] += pulse


        ## Calculate Ex from Dx
        for k in 1:KE
            ex[k] = ga[k] * (dx[k] - ix[k] - sx[k])
            ix[k] += gb[k]*ex[k]
            sx[k] = del_exp*sx[k] + gc[k]*ex[k]
        end

        ## Calculate the Fourier transform of Ex
        for k in 1:KE, m in 1:Nfreq
            real_pt[m, k] += cos(arg[m]*T) * ex[k]
            imag_pt[m, k] -= sin(arg[m]*T) * ex[k]            
        end

        ## Fourier transform of input pulse
        if T < 100
            for i in 1:Nfreq
                real_in[i] += cos(arg[i]*T) * ex[kc+5]
                imag_in[i] -= sin(arg[i]*T) * ex[kc+5]
            end
        end

        ## Boundary conditions
        ex[1] = ex_low_m2
        ex_low_m2 = ex_low_m1
        ex_low_m1 = ex[2]
        ex[KE] = ex_high_m2
        ex_high_m2 = ex_high_m1
        ex_high_m1 = ex[KE-1]

        ## Calculate the Hy field
        for k in 1:KE-1
            hy[k] += 0.5(ex[k] - ex[k+1])
        end

    end

    ## Write the E field out to a file
    open("Ex", "w") do f
        for k in 1:KE
            @printf(f, "  %6.3f \n", ex[k])
        end
    end

    ## Calculate the amplitude and phase of each frequency
    # Amplitude and phase of the input pulse
    for m in 1:Nfreq
        amp_in[m] = hypot(imag_in[m], real_in[m])
        phase_in[m] = atan2(imag_in[m], real_in[m])
        @printf("%d Input pulse: %8.4f %8.4f %8.4f %7.2f\n",
            m, real_in[m], imag_in[m], amp_in[m], (180.0/pi)*phase_in[m])
        for k in 1:KE
            ampn[m, k] = (1.0/amp_in[m])*hypot(real_pt[m, k], imag_pt[m, k])
            phasen[m, k] = atan2(imag_pt[m, k], real_pt[m, k]) - phase_in[m]
        end
    end

        #=
        for k in 1:KE
            @printf("%d %6.3f %6.3f %6.3f \n",
            k, ampn[1, k], ampn[2, k], ampn[3, k])
        end
        =#

    ## Write amplitude field out to a file "Amp"
    for m in 1:Nfreq
        open("Amp" * "$m", "w") do f
            for k in ampn[m, :]
                @printf(f, "    %8.5f \n", k)
            end
        end
    end
    

    ## Plot Ex and Hy
    figure(figsize=(15,3))
    plot([1:KE], ex)
    title("Ex")
    # figure()
    # plot([1:KE], hy)
    # title("Hy")

    ## Plot the Fourier amplitude
    for i in 1:Nfreq
        figure(figsize=(15,3))
        plot([1:KE], ampn[i,:]')
        title("Amp"*"$i")
    end
end


