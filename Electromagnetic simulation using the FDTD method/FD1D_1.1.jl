#= 
    1D FDTD simulation in free space
    ex[k] = ex[k] + 0.5 * (hy[k-1] - hy[k])
    hy[k] = hy[k] + 0.5 * (ex[k] - ex[k+1])
=#

using PyPlot

const KE = 200 # KE is the number of cells to be used
ex = zeros(KE)
hy = zeros(KE)

kc = KE /2     # Center of the problem space
t0 = 40.0      # Center of the incident pulse
spread = 12    # Width of the incident pulse
T = 0
NSTEPS = 1

while NSTEPS >= 1
    @printf("NSTEPS --> ")      # NSTEPS is the number of times the
    NSTEPS = int(readline())    # main loop has excuted
    @printf("%d \n", NSTEPS)
 
    ### Main FDTD loop
    for n = 1:NSTEPS
        T += 1                  # T keeps track of the total number
                                # of times the main loop is excuted
        
        # Calculate the Ex field
        for k = 2:KE
            ex[k] += 0.5(hy[k-1] - hy[k])
        end
        
        # Put a Gaussia pulse in the middle
        pulse = exp(-0.5((t0-T)/spread)^2.0)
        ex[kc] = pulse
        @printf("%5.1f  %6.2f\n", t0-T, ex[kc])
        
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
    open("Ew", "w") do f
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
    figure()
    plot([1:KE], ex)
    title("Ex")
    figure()
    plot([1:KE], hy)
    title("Hy")
end
