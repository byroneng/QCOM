program qcom

!     pgf90 -c print.f
!     pgf90 -o qcom QCOM.f90 print.o
!     ./qcom

!     v. 2 (Sept 2000) --
!     Array subscript ranges now start with 0 instead of 1.	

!     v. 3 (Feb 2009) --
!     All comments now start with !
!     Added f90 declarations and procedures.

      integer dj 
      real dk, ekth

      parameter (jt = 1, kt = 20)
      parameter (jv = jt, kv = kt)
      real, dimension (0:jv+1, 0:kv+1) :: v
      real, dimension (1:jv, 1:kv, 2) :: fv

!     Assign parameters for run

!      parameter (jw = jt, kw = kt)
      parameter (jth = jt, kth = kt)
!      parameter (jp = jt, kp = kt)

!     Initialize all values of arrays for v, w, theta, pi, fv, fw, ftheta, fpi,
!     including values for boundary points for v, w, theta, pi.

!      real, dimension (0:jw+1, 0:kw+1) :: w
!      real, dimension (1:jw, 1:kw, 2) :: fw

      real, dimension (0:jth+1, 0:kth+1) :: theta
      real, dimension (1:jth, 1:kth, 2) :: ftheta

!      real, dimension (0:jp+1, 0:kp+1) :: pi
!      real, dimension (1:jp, 1:kp, 2) :: fpi

      parameter (tmax = 200., dt = 1.25, Nout = 10)
      parameter (ITTMAX = tmax/dt)
!      real, dimension (1:FLOOR((ITTMAX * dt)/Nout), 0:kv+1) :: Vout
      real, dimension (1:FLOOR((ITTMAX * dt)/Nout), 0:kth+1) :: THout

      CALL INIT  ! initialize all variables

      open(98, file='qout.dat') ! Open the file 'qout.dat'

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

!      write(*,*) 'ITT = ',ITT

      CALL STEP ( N1, N2, A, B ) ! do subsequent time steps

      end do
	
!     END-OF-RUN OUTPUT ROUTINES GO HERE

 !     write(*,*) 'Calling PRINT'

!      CALL PRINT(Vout, FLOOR((ITTMAX * dt)/Nout), 1, ip, 1, kv+2, 'Vout')
      CALL PRINT(THout, FLOOR((ITTMAX * dt)/Nout), 1, ip, 1, kth+2, 'THout')

      close(98)   ! Close 'qout.dat'

!      write(*,*) 'J,1', theta(:,1)
!      write(*,*) '1,K', theta(1,:)

!      write(*,*) 'End of script, saved qout.dat'

contains

      SUBROUTINE STEP ( N1, N2, A, B )

 !     write(*,*) 'Called STEP'

!     This is the entire subroutine.

      CALL RCALC ( N2 ) ! calculate forcing terms from variables at current time
      CALL AB ( N1, N2, A, B ) ! update variables using a time scheme
      CALL BOUND ! apply boundary conditions to variables

      if (mod(ITT*DT,Nout*1.) .eq. 0) then
            ip = ip+1
!            Vout(ip,:) = v(1,:)
            THout(ip,:) = theta(1,:)
            write(98,*) theta(1,:)

      end if

!      write(*,*) theta
!
      END SUBROUTINE STEP

      SUBROUTINE RCALC (  N2 )

!      write(*,*) 'Called RCALC'

!     CALCULATES FORCING TERMS FOR V(J,K), ETC.; STORES THEM  IN FV(J,K,N2), ETC.

      DO K = 1, KT
      DO J = 1, JT
!      FV(J,K,N2)  = [fv for v(j,k)]
      ftheta(J,K,N2) = ekth*(((theta(J,k+1) - (2*theta(J,k))) + theta(J,k-1)) / (dk**2))
 !     write(*,*) 'k=',k,'ftheta = ', ((theta(J,k+1) - (2*theta(J,k))) + theta(J,k-1))
      END DO
      END DO
      
!      if (mod(ITT,Nout) .eq. 0) then
 !     if (ITT >= 140) then
 !           write(*,*) ITT,'k=1, j=1: ',theta(1,2),' - 2*',theta(1,1),' + ',theta(1,0)
 !           write(*,*) ' = ',((theta(1,2)-(2*theta(1,1)))+theta(1,0))/(dk**2)
 !     end if

!     ETC (forcing for w, theta, and pi)

      END SUBROUTINE RCALC

      SUBROUTINE AB ( N1, N2, A, B )

 !     write(*,*) 'Called AB'

!     THE FOLLOWING LOOP UPDATES V USING EITHER THE FORWARD OR THE ADAMS-BASHFORTH 
!     SCHEME DEPENDING ON THE VALUES OF A, B.  
!     SUBSCRIPT N2 OF FV ALWAYS REFERS TO THE MOST RECENTLY CALCULATED VALUES FOR FV.

      DO K = 1, KT
      DO J = 1, JT
!      V(J,K)  = V(J,K)  + DT * ( A * FV(J,K,N2)  + B * FV(J,K,N1) )
      theta(J,K) = theta(J,K) + (DT * (A*ftheta(J,K,N2) + B*ftheta(J,K,N1)))
!      write(*,*) 'AB Correction: ', (DT * (A*ftheta(J,K,N2) + B*ftheta(J,K,N1)))
      END DO
      END DO

!      write(*,*) 'ftheta(1,3) = ', ftheta(1,3,N2),'and',ftheta(1,3,N1), 'theta(1,3) = ', theta(1,3)

!     ETC (update w, theta, and pi)

      END SUBROUTINE AB

      SUBROUTINE BOUND

 !     write(*,*) 'Called BOUND'

!     apply boundary conditions for v

      do j = 1, jt
      theta(j,1) = 288. !Boundaries set to constant temperature
      theta(j,kt) = 288.-1.2 
      end do

      do j = 1, jt
      v(j,0) = v(j,1) ! free slip b.c.
      v(j,kt+1) = v(j,kt) ! free slip b.c.
      theta(j,0) = theta(j,1)
      theta(j,kt+1) = theta(j,kt)
      end do

      do k = 0, kt+1
      v(0,k) = v(jt,k) ! periodic b.c.
      v(jt+1,k) = v(1,k) ! periodic b.c.
      theta(0,k) = theta(jt,k)
      theta(jt+1,k) = theta(1,k)
      end do



!     ETC (apply b.c. for w, theta, and pi)

      END SUBROUTINE BOUND
      
      SUBROUTINE INIT
		
!      write(*,*) 'Called INIT'

!     initialize all variables 

      ekth = 100. !eddy viscosity
      dk = 25. !Vertical gridsize
      dj= 1 !gridsize (y-component)

      do k=1, kt
            theta(:,k) = (288. + (288. - (1.2)))/2.
      end do

      CALL BOUND

 !     write(*,*) 'Initial theta: ',theta(1,:)

      
      END SUBROUTINE INIT
      
end program qcom
