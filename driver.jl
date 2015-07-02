include("Circular Waveguide Mode Calculations/TE_ModeAnalysis.jl")

cond1 = ModeCondition(0, 3, 0.55, 95e9)

x, y, Emag = calTE(cond1)

Emag0 = vec(Emag[(1+end)/2, :])


include("find_localextrema.jl")
maxima, minima = find_localextrema(Emag0)
xmax = [x[i] for i in maxima]
