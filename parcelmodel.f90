program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ATMOS 6150
! Byron Eng
!
! This program performs a saturation adjustment on the parcel.
! Input: THETAstar Wstar Lstar Pnpp (before adjustment; after other processes)
! Output: THETAnpp Wnpp Lnpp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	implicit none

		real THETAstar, Wstar		! Initial Potential Temperature, Initial Water Vapor
		real     Lstar, Pnpp 		! Initial Liquid Water, Pressure (n+1)
		real        Po, Rv, Cp		! Ref Pressure, Gas Const, Spec Heat
		real      gamm, exnr		! Gamma, Exner Function
		real        Ws, es, alph	! Saturation Mixing Ratio, Saturation VP, alpha
		real         L, T, THETAnpp	! Latent Heat, Temperature, Potential Temp (n+1)
		real      Wnpp, Lnpp, R 	! Mixing Ratio (n+1), Liquid Water (n+1), Dry gas const
		real	THETAe, Wtot, RH	! Equivalent potential temp
		
		integer,dimension(5000)	:: Parray			! Pressure array
		character(15) filename
		real     dP, Pmin, Pstar	! Pressure increment, min pressure, pressure
		integer		 i

		print *, 'Enter maximum pressure (Pa)'
		read (*,*) Pstar

		print *, 'Enter minimum pressure (Pa)'
		read (*,*) Pmin

		print *, 'Enter initial temperature (K)'
		read (*,*) T

		print *, 'Enter initial water vapor concentration (g/g)'
		read (*,*) Wstar

		print *, 'Enter initial liquid water concentration (g/g)'
		read (*,*) Lstar

		dP    = (Pstar-Pmin)/size(Parray)					
		Po    = 100000
		Rv    = 461.5
		R     = 287
		Cp    = 1004
		L = 2.5e6
		Pnpp = Pstar - dP

		do i = 1,size(Parray)
			Parray(i) = Pstar - (dP*(i-1))
		end do

		filename = './parcelout.txt'
		write(*,*) 'Saving File as: '//filename
		open(98, file=filename)

		do i = 1,size(Parray)

			exnr = (Pnpp/Po)**(R/Cp)
			gamm  = L/(Cp*exnr)

			IF (i==1) THEN
			THETAstar    = T/exnr
			END IF

			T = THETAstar * exnr

			! Wexler's formula for es(T)
  			es = exp(   (-0.29912729e4  *((T)**(-2)))  + &
        	            (-0.60170128e4  *((T)**(-1)))  + &
                   		( 0.1887643854e2*((T)**( 0)))  + & 
                   		(-0.28354721e-1 *((T)**( 1)))  + &
                   		( 0.17838301e-4 *((T)**( 2)))  + &
                   		(-0.84150417e-9 *((T)**( 3)))  + &
                   		( 0.44412543e-12*((T)**( 4)))  + &
                   		( 0.2858487e1*log( T)))

			Ws = 0.622*(es/(Pnpp - es))

			alph = 0.622*((exnr*Pnpp)/(Pnpp - es)**2)*((L*es)/(Rv*(T**2)))
			
			!Saturated case
			THETAnpp = THETAstar + (gamm/(1+gamm*alph))*(Wstar - Ws)
			Wnpp = Ws + alph*(THETAnpp - THETAstar)
			Lnpp = Wstar + Lstar - Wnpp

			IF (dP > 0) THEN
			Lnpp = Lnpp - (.0002*dP*Lnpp) ! Precip. -dl/dp = -Cl  where C = .02 mb^-1
			END IF

			!Unsaturated case
			IF (Lnpp < 0) THEN
				Lnpp = 0
				Wnpp = Wstar + Lstar
				THETAnpp = THETAstar - gamm*(Wnpp - Wstar)
			END IF

			THETAstar = THETAnpp
			THETAe = THETAstar*exp((L*Ws)/(Cp*T))
			Wstar = Wnpp
			Lstar = Lnpp
			Pnpp = Pnpp - dP
			Wtot = Wstar + Lstar
			RH   = (Wstar/Ws)*100
			write(98,*) Pnpp,T,THETAstar,THETAe,Lstar,Wstar,Wtot,RH

			end do
		close(98)

end program