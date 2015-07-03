#=
    calculate rectangular TE10 field distribution at center
    and find locations of couplings
    Date:2015-07-03
    Reference: Pozar. Microwave Engineering.
=#
using PyPlot
f = 95e9 # frequency [Hz]
omega = 2pi*f # angular frequency
μ = 4e-7pi # permeability
a = 2.032e-3 # length of width[m]

N = 201 # number of points
x = linspace(0,a,N)
E(x) = -im*omega*μ*a/pi*sin(pi*x/a) # TE10 mode distribution


ratio = [1, 1.2659, 2.1292] # ratio of TE03 maxima
total = abs(quadgk(E, x[1], x[end])[1])
areas = zeros(N)
for i in 1:N
    areas[i] = abs(quadgk(E, x[i], x[end])[1])
end

proportion = areas/total

firstmaxratio = ratio[3]/sum(ratio)
secondmaxratio = ratio[2]/sum(ratio)

firstmaxlocind = indmin(abs(proportion - firstmaxratio))

# location where the first coupling is done
firstmaxloc = x[firstmaxlocind]

areas2 = zeros(firstmaxlocind)
for i in 1:firstmaxlocind
    areas2[i] = (abs(quadgk(E, x[i], x[firstmaxlocind])[1]))
end

proportion = areas2/total
secondmaxlocind = indmin(abs(proportion - secondmaxratio))
secondmaxloc = x[secondmaxlocind]
