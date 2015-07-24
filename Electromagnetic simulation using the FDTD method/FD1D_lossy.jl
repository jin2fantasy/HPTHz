#= 
    1D FDTD simulation in lossy dielectric
    main equations:
    ex[k] = ca[k]*ex[k] + cb[k]*(hy[k-1] - hy[k])
    hy[k] = hy[k] + 0.5*(ex[k] - ex[k+1])
    where eaf = dt*σ/(2ε0*εr)
    ca[k] = (1.0 - eaf)/(1.0 + eaf)
    cb[k] = 0.5/(εr*(1.0 + eaf))
=#

using PyPlot
const c = 299792458.0  # speed of light [m/s]
const μ0 = 4e-7pi      # vacuum permeability
const ε0 = 1/(μ0*c^2) # vacuum permittivity
const εr = 4.0        # relative permittivity of dielectric
const σ = 0.04         # conductivity of dielectric
const KE = 200         # KE is the number of cells to be used
ex = zeros(KE)
hy = zeros(KE)
ca = ones(KE)         # parameter 1
cb = fill(0.5, KE)    # parameter 2
t0 = 40.0             # Center of the incident pulse (gaussian pulse)
spread = 12           # Width of the incident pulse  

freq_in = 700e6       # frequency [Hz]
λ = c/freq_in        # wavelength [m]
Nsampling = 10        # number of points per wavelength
ddx = λ/Nsampling    # cell size
dt = ddx/2c           # Calculate the time step

eaf = dt*σ/(2ε0*εr)
dielectric_range = 0.5KE:KE # range of dielctric
ca[dielectric_range] = (1.0 - eaf)/(1.0 + eaf) # conductivity
cb[dielectric_range] = 0.5/(εr*(1.0 + eaf))   # dielectric constants

T = 0
NSTEPS = 1
# parameters for absorbing boundary condition
ex_low_m1 = ex_low_m2 = ex_high_m1 = ex_high_m2 = 0.0

while NSTEPS >= 1
    @printf("NSTEPS --> ")      # NSTEPS is the number of times the
    NSTEPS = int(readline())    # main loop has excuted
    @printf("%d \n", NSTEPS)
    NSTEPS >= 1 || (break)

    ### Main FDTD loop
    for n = 1:NSTEPS
        T += 1                  # T keeps track of the total number
                                # of times the main loop is excuted
        
        # Calculate the Ex field
        for k = 2:KE
            ex[k] = ca[k]*ex[k] + cb[k]*(hy[k-1] - hy[k])
        end
        
        # Gaussia pulse
        # pulse = exp(-0.5((t0-T)/spread)^2.0)
        # Sinusoidal pulse
        pulse = sin(2pi*freq_in*dt*T)
        ex[5] += pulse
        
        # Absorbing boundary condition (valid for free space)
        ex[1]      = ex_low_m2
        ex_low_m2  = ex_low_m1
        ex_low_m1  = ex[2]
        ex[KE]     = ex_high_m2
        ex_high_m2 = ex_high_m1
        ex_high_m1 = ex[KE-1]

        # Calculate the Hy field
        for k = 1:KE-1
            hy[k] += 0.5(ex[k] - ex[k+1])
        end
        
    end  ### End of Main FDTD loop
    
    # At the end of the calculation, print out
    # the Ex and Hy fields
    for k = 1:KE
        @printf("%3d   %6.2f   %6.2f\n", k, ex[k], hy[k])
    end
    
    # Write the E field out to a file "Ex"
    open("Ex", "w") do f
        for k = 1:KE
            @printf(f, "  %6.2f \n", ex[k])
        end
    end
    
    # Write the H field out to a file "Hy"
    open("Hy", "w") do f
        for k = 1:KE
            @printf(f, "  %6.2f \n", hy[k])
        end
    end
    
    @printf("T = %5.0f\n", T)
    
    # Plot Ex and Hy
    figure(figsize=(15,3))
    plot([1:KE], ex)
    title("Ex")
    # figure()
    # plot([1:KE], hy)
    # title("Hy")
end
