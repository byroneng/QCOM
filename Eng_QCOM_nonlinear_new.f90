program qcom

!     pgf90 -c print.f
!     pgf90 -o qcom QCOM.f90 print.o
!     ./qcom

!     v. 2 (Sept 2000) --
!     Array subscript ranges now start with 0 instead of 1.	

!     v. 3 (Feb 2009) --
!     All comments now start with !
!     Added f90 declarations and procedures.

      real dk, ekth, ekp, ekv, g
      real Cp, th0, Cs, H, L, dj

      parameter (jt = 20, kt = 10)
      parameter (jv = jt, kv = kt)
      real, dimension (0:jv+1, 0:kv+1) :: v
      real, dimension (1:jv, 1:kv, 2) :: fv

!     Assign parameters for run

      parameter (jw = jt, kw = kt)
      parameter (jth = jt, kth = kt)
      parameter (jp = jt, kp = kt)

!     Initialize all values of arrays for v, w, theta, pi, fv, fw, ftheta, fpi,
!     including values for boundary points for v, w, theta, pi.

      real, dimension (0:jw+1, 0:kw+1) :: w
      real, dimension (1:jw, 1:kw, 2) :: fw

      real, dimension (0:jth+1, 0:kth+1) :: thetao
      real, dimension (0:jth+1, 0:kth+1) :: theta
      real, dimension (1:jth, 1:kth, 2) :: ftheta

      real, dimension (0:jp+1, 0:kp+1) :: pi
      real, dimension (1:jp, 1:kp, 2) :: fpi

      parameter (tmax = 100000., dt = .1, Nout = 10)
      parameter (ITTMAX = tmax/dt)
      real, dimension (1:FLOOR(tmax/Nout), 0:kv+1) :: Vout
      real, dimension (1:FLOOR(tmax/Nout), 0:kth+1) :: THout

      CALL INIT  ! initialize all variables

      ITT = 1 ! itt is time step index

!     USE FORWARD SCHEME TO do first step

      A = 1.
      B = 0.
      N1 = MOD ( ITT    , 2 ) + 1
      N2 = MOD ( ITT - 1, 2 ) + 1

      CALL STEP ( N1, N2, A, B ) ! do first time step

!     ADAMS - BASHFORTH TWO - LEVEL SCHEME

      A =   3. / 2. 
      B = - 1. / 2. 

      ITTNOW = ITT  + 1

      DO ITT = ITTNOW, ITTMAX

      N1 = MOD ( ITT    , 2 ) + 1   !N1 and N2 alternnate values of 1 and 2 here
      N2 = MOD ( ITT - 1, 2 ) + 1


      CALL STEP ( N1, N2, A, B ) ! do subsequent time steps

      end do
	
!     END-OF-RUN OUTPUT ROUTINES GO HERE

            open(91, file='v.dat') ! Open the file 'v.dat'
                  do k=0, kt+1
                        write(91,*) v(:,k)
                  end do
            close(91)


            open(92, file='w.dat') ! Open the file 'w.dat'
                  do k=0, kt+1
                        write(92,*) w(:,k)
                  end do
            close(92)


            open(93, file='theta.dat') ! Open the file 'theta.dat'
                  do k=0, kt+1
                        write(93,*) theta(:,k)-thetao(:,k)
                  end do
            close(93)


           open(94, file='pi.dat') ! Open the file 'pi.dat'
                  do k=0, kt+1
                        write(94,*) pi(:,k)
                  end do
           close(94)

      write(*,*) 'End of script, saved v, w, theta, and pi .dat files'



contains

      SUBROUTINE STEP ( N1, N2, A, B )


!     This is the entire subroutine.

      CALL RCALC ( N2 ) ! calculate forcing terms from variables at current time
      CALL AB ( N1, N2, A, B ) ! update variables using a time scheme
      CALL BOUND ! apply boundary conditions to variables

!      if (mod(ITT*DT,Nout*1.) .eq. 0) then
!            ip = ip+1
!            Vout(ip,:) = v(10,:)
!            THout(ip,:) = theta(10,:)
!            write(98,*) v(1,:), w(1,:), theta(1,:), pi(1,:)
!      end if

!!!!!!!!!!!!!!!!!!!!!!!! Print step 5 (for debugging) !!!!!!!!!!!!!!!!!!!!!!!!!
!       if (itt == 5) then
!             open(91, file='v.dat') ! Open the file 'v.dat'
!                   do k=1, kt
!                         write(91,*) v(:,k)
!                   end do
!             close(91)


!             open(92, file='w.dat') ! Open the file 'w.dat'
!                   do k=1, kt
!                         write(92,*) w(:,k)
!                   end do
!             close(92)


!             open(93, file='theta.dat') ! Open the file 'theta.dat'
!                   do k=1, kt
!                         write(93,*) theta(:,k)
!                         write(*,*) theta(:,k)
!                   end do
!             close(93)


!            open(94, file='pi.dat') ! Open the file 'pi.dat'
!                   do k=1, kt
!                         write(94,*) pi(:,k)
!                   end do
!            close(94)
!       end if
!!!!!!!!!!!!!!!!!!!!!!!! Print step 5 (for debugging) !!!!!!!!!!!!!!!!!!!!!!!!!

!       if (itt < 100) then
!             write(*,*) 'fv =', fv(5,5,2), &
!                   'fw =', fw(5,5,2), &
!                   ' ftheta = ', ftheta(5,5,2), &
!                   'fpi =', fpi(5,5,2)

!       end if

      if (itt<100) then
!            write(*,*) 'v vert = ', v(5,:)
!            write(*,*) 'w vert = ', w(5,:)
            write(*,*) 'itt=',itt,'pi vert= ', pi(5,:)

      end if

      END SUBROUTINE STEP

      SUBROUTINE RCALC (  N2 )


!     CALCULATES FORCING TERMS FOR V(J,K), ETC.; STORES THEM  IN FV(J,K,N2), ETC.

      DO K = 1, KT
      DO J = 1, JT
      fv(J,K,N2)  = - (v(j,k)*((v(j+1,k)-v(j-1,k))/(2*dj)))                    &
                    - (((((w(j,k)+w(j+1,k))/2)*((v(j,k+1)-v(j,k))/dk))        &
                     +(((w(j,k-1)+w(j+1,k-1))/2)*((v(j,k)-v(j,k-1))/dk))/2)) &
                    - (Cp*thetao(j,k)*(pi(j+1,k) - pi(j,k))/dj)                &
                    + (ekv*(v(j,k+1) - (2*v(j,k)) + v(j,k-1)) / (dk**2))      &
                    + (ekv*(v(j+1,k) - (2*v(j,k)) + v(j-1,k)) / (dj**2))
      END DO
      END DO


      DO K = 1, KT
      DO J = 1, JT
      fw(j,k,N2) = -(((((v(j,k+1)+v(j,k))/2)*((w(j+1,k)-w(j,k))/dj))               &
                        +(((v(j-1,k+1)+v(j-1,k))/2)*((w(j,k)-w(j-1,k))/dj)))/2)    &
            - (w(j,k)*((w(j,k+1)-w(j,k-1))/(2*dk)))                                  &
            - (Cp*(0.5*(thetao(j,k+1)+thetao(j,k)))*(pi(j,k+1) - pi(j,k))/dk)        &
            + (g*((((theta(j,k+1) + theta(j,k))/2)/((thetao(j,k+1)+thetao(j,k))/2))-1))  &
            + (ekp*(w(j,k+1) - (2*w(j,k)) + w(j,k-1)) / (dk**2))                    &
            + (ekp*(w(j+1,k) - (2*w(j,k)) + w(j-1,k)) / (dj**2))
      END DO
      END DO


      DO K = 1, KT
      DO J = 1, JT
      ftheta(J,K,N2) =   - (((v(j,k)*((theta(j+1,k)-theta(j,k))/dj))                          &
                           +(v(j-1,k)*((theta(j,k)-theta(j-1,k))/dj)))/2)                   &
                        + (ekth*(theta(J,k+1) - (2*theta(J,k)) + theta(J,k-1)) / (dk**2))   &
                        + (ekth*(theta(j+1,k) - (2*theta(j,k)) + theta(j-1,k)) / (dj**2))   &
                        - (0.5*(w(j,k)  *(theta(j,k+1) - theta(j,k  ))                     &
                            + (w(j,k-1)*(theta(j,k  ) - theta(j,k-1)))) / dk)
      END DO
      END DO


      DO K = 1, KT
      DO J = 1, JT
      fpi(j,k,N2) = -((Cs**2)/(Cp*(thetao(j,k)**2)))                              &
                  * ((((thetao(j,k)*v(j,k)) - (thetao(j,k)*v(j-1,k))) / dj)      &
                  + (((thetao(j,k+1) + thetao(j,k))*w(j,k)/2) - ((thetao(j,k) + thetao(j,k-1))*w(j,k-1)/2))/dk)

      END DO
      END DO


!     ETC (forcing for w, theta, and pi)

      END SUBROUTINE RCALC

      SUBROUTINE AB ( N1, N2, A, B )


!     THE FOLLOWING LOOP UPDATES V USING EITHER THE FORWARD OR THE ADAMS-BASHFORTH 
!     SCHEME DEPENDING ON THE VALUES OF A, B.  
!     SUBSCRIPT N2 OF FV ALWAYS REFERS TO THE MOST RECENTLY CALCULATED VALUES FOR FV.

      DO K = 1, KT
      DO J = 1, JT
      V(J,K)     = V(J,K)     + (DT  * (A * FV(J,K,N2)     + B * FV(J,K,N1)))
      w(j,k)     = w(j,k)     + (DT  * (A * fw(j,k,N2)     + B * fw(j,k,N1)))
      theta(J,K) = theta(J,K) + (DT  * (A * ftheta(J,K,N2) + B * ftheta(J,K,N1)))
      pi(j,k)    = pi(j,k)    + (DT  * (A * fpi(j,k,N2)    + B * fpi(j,k,N1)))
      END DO
      END DO


!     ETC (update w, theta, and pi)

      END SUBROUTINE AB

      SUBROUTINE BOUND


!     apply boundary conditions for v

      do j = 1, jt 
      v(j,0) = v(j,1) ! free slip b.c.
      v(j,kt+1) = v(j,kt) ! free slip b.c.
      w(j,0) = 0 ! no vert. motion at surface
      w(j,kt+1) = 0.
      theta(j,0) = thetao(j,0)
      theta(j,kt+1) = thetao(j,kt+1)
      pi(j,0) = pi(j,1)
      pi(j,kt+1) = pi(j,kt)
      end do

      do k = 0, kt+1
      v(0,k) = v(jt,k) ! periodic b.c.
      v(jt+1,k) = v(1,k) ! periodic b.c.
      w(0,k) = w(jt,k)
      w(jt+1,k) = w(1,k)
      theta(0,k) = theta(jt,k)
      theta(jt+1,k) = theta(1,k)
      pi(0,k) = pi(jt,k)
      pi(jt+1,k) = pi(1,k)
      end do



!     ETC (apply b.c. for w, theta, and pi)

      END SUBROUTINE BOUND
      
      SUBROUTINE INIT
		
!     initialize all variables 

      ekth = 50. !eddy viscosity
      ekv  = 50.
      ekp  = 50.
      H = 500.
      L = (2**(3/2))*H
      delth = 1.2
      Cs = 50.
      dk = H/kt !Vertical gridsize
      dj = L/jt !y- gridsize
      Cp = 1004.
      g = 9.81


      do k=0, kt+1
            thetao(:,k) = (288. - (delth*k*dk/H))
            theta(:,k) = (288. - (delth*k*dk/H))
            v(:,k) = 0.0
            w(:,k) = 0.0
            pi(:,k) = 0.0
      end do


      do j=1, jt
            theta(j,floor(kt/2.):ceiling((kt+.5)/2.)) = theta(j,floor(kt/2.):ceiling((kt+.5)/2.)) + &
                                                            (.05*2.4)*cos((2*3.14159/L)*j*dj)
      end do
    
      CALL BOUND
      
      END SUBROUTINE INIT
      
end program qcom
