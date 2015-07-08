# Original script from Matlab TE_ModeAnalysis.m

# The cylindrical waveguide electric field
# The code calculates linear combination of two linearly polarized modes
# with arbitrary phase. For the ideal circulary polarized mode, the field
# pattern comprises of a linear combination of two linearly polarized modes
# with pi/2 phase difference. Nov. 2013, EMC

using Roots
using Calculus
using PyPlot

immutable ModeCondition
	m::Int # azimuthal number
	n::Int # radial number
	Rw::Float64 # [cm] radius of the cavity
	f::Float64 # [Hz] frequency
end

function BesselJPrimeRoots(m::Int, n::Int)
    limit = 50
    besseljm(x) = besselj(m, x)
    dbesseljm = derivative(besseljm)
    dbjzeros = fzeros(dbesseljm, 0.0, limit) # limit can be changed
                                       # if n is high
    if dbjzeros[1] == 0.0
        deleteat!(dbjzeros, 1)
    end

    testonemore(x::Float64) = dbesseljm(x)*dbesseljm(nextfloat(x)) > 0 &&
                                dbesseljm(x)*dbesseljm(prevfloat(x)) > 0

    while testonemore(dbjzeros[1])
        deleteat!(dbjzeros, 1)
    end
    n > length(dbjzeros) && (error("n is too high; increase limit"))
    dbjzeros[n]
end


function calTE(cond::ModeCondition)
    m = cond.m
    n = cond.n
    Rw = cond.Rw
    Rlim  = 1.2 * Rw # [cm] y axes limit in plot
    c = 2.997e10 # [cm/s] speed of light
    omega = cond.f * 2pi
	N = 213 # number of grid
	A = 1 	# amplitude of linearly pol mode
	# A = 1.8 # amplitude of conuter-(or co-) rotating modes
	B = 1	# Amplitude for linearly pol mode 1 (set it 1)
	phi1 = pi/2 #= phase difference, it should be pi/2 for ideal
				   circular rotation case =#
	phi2 = 0 	# angle of phase diff (set it 0)
	phi3 = 0 	# -10*pi/180 # angle of field rotation, for rotation matrix

	nump = BesselJPrimeRoots(m, n)
	kr = nump/Rw
	k = omega/c
	kz = sqrt(k^2 - kr^2)

	x = linspace(-Rlim,Rlim,N)
	y = linspace(-Rlim,Rlim,N)

    Ex = zeros(Complex, N, N)
	Ey = zeros(Complex, N, N)
	Exold = Eyold = 0im
    Exx = zeros(Complex, N)
	Eyy = zeros(Complex, N)

	for kk = 1:N
        for qq = 1:N
    		r = sqrt(x[kk]^2 + y[qq]^2)
    		phi = angle(complex(x[kk], y[qq]))
            temp1 = nump*r/Rw
    		Jmprime = (besselj(m-1, temp1) - besselj(m+1, temp1)) / 2
            Jm = besselj(m, temp1)

            if r <= Rw
                # from above equations, add arbitrary amplitude (A)
                # & phase (phi1)
                Exold = ((im*m/kr/r)*B*exp(im*phi2)*cos(m*phi)*cos(phi)*Jm
				 	+ im*B*exp(im*phi2)*sin(m*phi)*sin(phi)*Jmprime
                    - A*exp(im*phi1)*(im*m/kr/r)*sin(m*phi)*cos(phi)*Jm
                    + im*A*exp(im*phi1)*cos(m*phi)*sin(phi)*Jmprime)

                Eyold = ((im*m/kr/r)*B*exp(im*phi2)*cos(m*phi)*sin(phi)*Jm
				 	- im*B*exp(im*phi2)*sin(m*phi)*cos(phi)*Jmprime
                    - A*exp(im*phi1)*(im*m/kr/r)*sin(m*phi)*sin(phi)*Jm
                    - im*A*exp(im*phi1)*cos(m*phi)*cos(phi)*Jmprime)
            else
                Exold = complex(0.01)
                Eyold = complex(0.01)
            end
            Exx[qq] = Exold
            Eyy[qq] = Eyold
        end
        Ex[:, kk] = Exx[:]
        Ey[:, kk] = Eyy[:]
    end

    Ex1 = Ex.*complex(cos(phi3)) - complex(sin(phi3)).*Ey # rotation of field
														  # w.r.t.
                                        				  # observation plane
    Ey1 = complex(sin(phi3)).*Ex + Ey.*complex(cos(phi3))
    Ex = Ex1
    Ey = Ey1

    # Efield = Ex + exp(im*phi3)*Ey
    Exdb = 20 .* log(float64(abs(Ex)./maximum(abs(Ex))))
    Eydb = 20 .* log(float64(abs(Ey)./maximum(abs(Ey))))

    Emag = sqrt(float64((abs2(Ex)) + (abs2(Ey))))
    Edb = 20 .* log(Emag./maximum(Emag))

    Exphase = angle(Ex)
    Eyphase = angle(Ey)

    figure(1)
    subplot(221)
    img1 = imshow(abs(Ex), extent=[x[1],x[end],y[1],y[end]]);
    grid("on"); colorbar(img1); axis("on")
    subplot(222)
    img2 = imshow(abs(Ey), extent=[x[1],x[end],y[1],y[end]]);
    grid("on"); colorbar(img2); axis("on")
    subplot(223)
    img3 = imshow(Exphase, extent=[x[1],x[end],y[1],y[end]]);
    grid("on"); colorbar(img3); axis("on")
    subplot(224)
    img4 = imshow(Eyphase, extent=[x[1],x[end],y[1],y[end]]);
    grid("on"); colorbar(img4); axis("on")

    figure(2)
    subplot(221)
    clims = [-40 0];
    img5 = imshow(Exdb,vmin=-40,vmax=0, extent=[x[1],x[end],y[1],y[end]]);
    grid("on"); colorbar(img5); axis("on")
    subplot(222)
    img6 = imshow(Eydb,vmin=-40,vmax=0, extent=[x[1],x[end],y[1],y[end]]);
    grid("on"); colorbar(img6); axis("on")
    subplot(223)
    img7 = imshow(Edb,vmin=-40,vmax=0, extent=[x[1],x[end],y[1],y[end]]);
    grid("on"); colorbar(img7); axis("on")
    subplot(224)
    img8 = imshow(Emag, extent=[x[1],x[end],y[1],y[end]]);
    grid("on"); colorbar(img8); axis("on")

	return x, y, Emag
end
