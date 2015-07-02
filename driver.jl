include("Circular Waveguide Mode Calculations/TE_ModeAnalysis.jl")

cond1 = ModeCondition(0, 3, 0.55, 95e9)

calTE(cond1)
