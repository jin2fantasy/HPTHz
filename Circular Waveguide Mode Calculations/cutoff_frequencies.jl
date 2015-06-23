## Calculate cutoff frequencies of modes in the circular waveguide
## Date:2015-06-22
## Author:wonjin
## Reference: David. M. Pozar. Microwave Enginnering. 4th edition. p.124

Diam = 11e-3 # diameter of circular waveguide, in [m]
const c = 299792458 # speed of light, in [m/s]

using Calculus
# functions for derivatives of bessel functions
#=! a variable x shouldn't be declared before !=#
function deriv_besselj0(n::Number)
    dbj0 = differentiate("besselj0(x)", :x)
    global x = n
    return eval(dbj0)
end
function deriv_besselj1(n::Number)
    dbj1 = differentiate("besselj1(x)", :x)
    global x = n
    return eval(dbj1)
end

using Roots
# find the roots of derivatives of bessel functions
dbj0roots = fzeros(deriv_besselj0,-1,20)
dbj1roots = fzeros(deriv_besselj1,-1,20)

# cutoff frequencies for TE11, TE01, TE02, TE03 modes, respectively, in [Hz]
fc11 = c*dbj1roots[1]/(pi*Diam)
fc01 = c*dbj0roots[2]/(pi*Diam)
fc02 = c*dbj0roots[3]/(pi*Diam)
fc03 = c*dbj0roots[4]/(pi*Diam)
