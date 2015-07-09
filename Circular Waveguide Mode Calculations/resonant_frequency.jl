using Calculus
using Roots

function BesselJPrimeRoots(m::Int, n::Int)
    besseljm(x) = besselj(m, x)
    dbesseljm = derivative(besseljm)
    dbjzeros = fzeros(dbesseljm, 0.0, 50) # limits can be changed
                                       # if n is high
    if dbjzeros[1] == 0.0
        deleteat!(dbjzeros, 1)
    end

    testonemore(x::Float64) = dbesseljm(x)*dbesseljm(nextfloat(x)) > 0 &&
                                dbesseljm(x)*dbesseljm(prevfloat(x)) > 0

    while testonemore(dbjzeros[1])
        deleteat!(dbjzeros, 1)
    end

    dbjzeros[n]
end


const c = 299792458
const f = 95e9
const lambda0 = c/f
const a = 6.55e-3
const d = 4*lambda0
function resonant_frequency(m::Int, n::Int, l::Int)
    μᵣ = 1
    ϵᵣ = 1
    c/(2pi * √(μᵣ*ϵᵣ)) * √((BesselJPrimeRoots(m, n)/a)^2 + (l*pi/d)^2)
end

open("resonant_result.txt", "w") do f
    for m in 0, n in 1:4, l in 0:10
        print(f, "f$(m)$(n)$(l) is ")
        print(f, resonant_frequency(m, n, l)/1.0e9)
        println(f, " GHz")
    end
end
