# Original script from Matlab TE_ModeAnalysis.m

# The cylindrical waveguide electric field
# The code calculates linear combination of two linearly polarized modes
# with arbitrary phase. For the ideal circulary polarized mode, the field
# pattern comprises of a linear combination of two linearly polarized modes
# with pi/2 phase difference. Nov. 2013, EMC

using Roots
using Calculus
using PyPlot

type ModeCondition
	m::Int # azimuthal number
	n::Int # radial number
	Rw # [cm] radius of the cavity
	f # [Hz] frequency
end

function BesselJPrimeRoots(m::Int, n::Int)
	besseljm(x) = besselj(m, x)
	dbesseljm = derivative(besseljm)
	dbjzeros = fzeros(dbesseljm, 0, 20) # limits can be changed
										# if n is high
	dbjzeros[n+1]
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

	nump = BesselJPrimeRoots(cond.m, cond.n)
	kr = nump/Rw
	k = omega/c
	kz = sqrt(k^2 - kr^2)
	
	x = linspace(-Rlim,Rlim,N)
	y = linspace(-Rlim,Rlim,N)

	Ex = Ey = Exold = Eyold = Exx = Eyy = zeros(Complex, N)

	for kk = 1:N
        for qq = 1:N
    		r = sqrt(x[kk].^2 + y[qq].^2)
    		phi = angle(complex(x[kk], y[qq]))
            temp1 = nump.*r./Rw
    		Jmprime = (besselj(m-1, temp1) - besselj(m+1, temp1)) ./ 2
            Jm = besselj(m, temp1)

            if r <= Rw
                # from above equations, add arbitrary amplitude (A)
                # & phase (phi1)
                Exold = ((im.*m./kr./r).*B*exp(im*phi2)*cos(m.*(phi)).*cos(phi)
                    .* Jm + im.*B*exp(im*phi2)*sin(m.*(phi)).*sin(phi).*Jmprime
                    - A.*exp(im*phi1).*(im.*(m)./kr./r).*sin((m).*(phi)).*cos(phi)
                    .*Jm + im.*A.*exp(im*phi1).*cos((m).*(phi)).*sin(phi)
                    .*Jmprime)
                Eyold = ((im.*m./kr./r).*B*exp(im*phi2)*cos(m.*(phi)).*sin(phi)
                    .* Jm - im.*B*exp(im*phi2)*sin(m.*(phi)).*cos(phi).*Jmprime
                    - A.*exp(im*  phi1).* (im.*(m)./kr./r).*sin((m).*(phi))
                    .*sin(phi).*Jm  - im.*A.*exp(im*phi1).*cos((m).*(phi))
                    .*cos(phi).*Jmprime)  
            else
                Exold=10^(-2.0);
                Eyold=10^(-2.0);
            end
            Exx[qq] = Exold
            Eyy[qq] = Eyold
        end
        Ex = [Ex Exx]
        Ey = [Ey Eyy]
    end

    Ex1 = Ex.*complex(cos(phi3)) - complex(sin(phi3)).*Ey # rotation of field w.r.t. 
                                        # observation plane
    Ey1 = complex(sin(phi3)).*Ex + Ey.*complex(cos(phi3))
    Ex = Ex1
    Ey = Ey1
    # Efield = Ex + exp(im*phi3)*Ey
    Exdb = 20 .* log(abs(Ex)./maximum(abs(Ex)))
    Eydb = 20 .* log(abs(Ey)./maximum(abs(Ey)))

    Emag = sqrt(((real(Ex)).^2 + (imag(Ex)).^2) + ((real(Ey).^2 + imag(Ey).^2)))
    Edb = 20 .* log(Emag./maximum(Emag))

    Exphase = angle(Ex)
    Exphase = angle(Ey)

    figure(1)
    subplot(221)
    imagesc(x,y,abs(Ex)); grid("on"); colorbar; axis("xy")
    subplot(222)
    imagesc(x,y,abs(Ey)); grid("on"); colorbar; axis("xy")
    subplot(223);
    imagesc(x,y,Exphase); grid("on"); colorbar; axis("xy")
    subplot(224);
    imagesc(x,y,Eyphase); grid("on"); colorbar; axis("xy")

    figure(2)
    subplot(221)
    clims = [-40 0];
    imagesc(x,y,Exdb,clims); grid("on"); colorbar; axis("xy")
    subplot(222)
    imagesc(x,y,Eydb,clims); grid("on"); colorbar; axis("xy")
    subplot(223)
    imagesc(x,y,Edb,clims); grid("on"); colorbar; axis("xy")
    subplot(224)
    imagesc(x,y,Emag); grid("on"); colorbar; axis("xy")

end