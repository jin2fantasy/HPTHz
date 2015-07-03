include("TE_ModeAnalysis.jl")

cond1 = ModeCondition(0, 3, 0.55, 95e9)

x, y, Emag = calTE(cond1)

Emag0 = vec(Emag[(1+end)/2, :])


include("find_localextrema.jl")
maxima, minima = find_localextrema(Emag0)
xmax = [x[i] for i in maxima]
Emax = [Emag0[i] for i in maxima]

println("maxima locations and their field values are")
for i = 1:length(xmax)
    @printf("% 0.5f   % 0.5f\n", xmax[i], Emax[i])
end
