using SymPy
a, x, b = symbols("a, x, b", real=true)
y = a*cos(0.4x - b) + a/2
dy = diff(y, x)
rotmatrix(θ) = [cosd(θ) -sind(θ); sind(θ) cosd(θ)]
rotatedvector = rotmatrix(45)*[0.382683, -0.923880]
solve([Eq(subs(dy, x, 0), (-0.382684/0.92388)), Eq(subs(y, x, 0), 2.592137)])

