program qcom

!     gfortran -c -fdefault-real-8 adjust.f90
!     gfortran -c -fdefault-real-8 wexler.f90
!     gfortran -c -fdefault-real-8 lowe.f90
!     gfortran -o qcom -fdefault-real-8 -Wall -fcheck=all Eng_QCOM.f90 adjust.o wexler.o lowe.o

!       real THSTAR, TH1, QVSTAR, QV1
!       parameter (HLF = 2500000., pzero = 100000.)
      implicit none

      real dk, ekth, ekp, ekv, La, esat
      real th0, Cs, H, L, dj, qvs, T, RELHUM


      !!! Parameters
      integer jt, kt, jv, kv, jw, kw, jth, kth, jp, kp, ittmax, Nout, j, k, ITTNOW
      integer itt, N1, N2
      real Cp, g, tmax, dt, delth, a, b, rgas
      real thetaadj, qvadj, qcadj, piadj, stabht, stabdpt, stabty


      logical animate, debug, cloudtxt
      parameter (jt = 40, kt = 100, jv = jt, kv = kt) 

      parameter (Cp=1005., g=9.8, rgas=287.04)
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
      real, dimension (0:jth+1, 0:kth+1) :: thetav
      real, dimension (0:jth+1, 0:kth+1) :: thetavo
      real, dimension (0:jth+1, 0:kth+1) :: thetal
      real, dimension (1:jth, 1:kth, 2) :: fthetal

      real, dimension (0:jp+1, 0:kp+1) :: pio
      real, dimension (0:jp+1, 0:kp+1) :: pi
      real, dimension (1:jp, 1:kp, 2) :: fpi

      real, dimension (0:jth+1, 0:kth+1) :: qv
      real, dimension (0:jth+1, 0:kth+1) :: qvo
      real, dimension (0:jth+1, 0:kth+1) :: qc
      real, dimension (0:jth+1, 0:kth+1) :: qw
      real, dimension (1:jth, 1:kth, 2) :: fqw

      parameter (tmax = 3000., dt = .1) 
      parameter (ITTMAX = int(tmax/dt), Nout = 10)

!!!!!!!!!!!!!!!!!!!! Initialize variables !!!!!!!!!!!!!!!!!!!!

      debug = .true.
      animate = .true.
      cloudtxt = .false.
      ekth = 50. !eddy viscosity
      ekv  = 50.
      ekp  = 50.
      H = 5000.
      L = (2.**(3./2.))*H/5.
      delth = .010 ! Lapse rate of bottom layer (K/m)
      stabty = .01 ! Lapse rate of stable layer (K/m)
      stabht = 1000. ! (un)Stable layer height (m)
      stabdpt = 500. ! (un)Stable layer depth (m)
      Cs = 50.
      dk = H/real(kt) !Vertical gridsize
      dj = L/real(jt) !y- gridsize
      La = 2.5e6 !J K^-1 kg^1
!      RELHUM = .99 !Initial relative humidity

      !Initial theta profile
      !Bottom layer
      do k= 0, floor(stabht/dk)
            thetao(:,k) = (288. - (delth*k*dk))
            qvo(:,k) = .005
      end do

      !Capping Inversion layer
      do k = floor(stabht/dk)+1, floor((stabht+stabdpt)/dk)
            thetao(:,k) = thetao(:,floor(stabht/dk)) + (stabty*((k*dk)-stabht))
            qvo(:,k) = .002
      end do

      !Neutral layer
      do k = floor((stabht+stabdpt)/dk), kt+1
            thetao(:,k) = thetao(:,floor((stabht+stabdpt)/dk))
            qvo(:,k) = 0.0001
      end do


      do k=0, kt+1
            theta(:,k) = thetao(:,k)!(288. - (delth*k*dk))
            v(:,k) = 0.0
            w(:,k) = 0.0
            pi(:,k) = 0.0
            qc(:,k) = 0.0
      end do

      do k=1, kt
            ftheta(:,k,1) = 0.0
            fthetal(:,k,1) = 0.0
            fv(:,k,1) = 0.0
            fw(:,k,1) = 0.0
            fpi(:,k,1) = 0.0

            ftheta(:,k,2) = 0.0
            fthetal(:,k,2) = 0.0
            fv(:,k,2) = 0.0
            fw(:,k,2) = 0.0
            fpi(:,k,2) = 0.0
      end do

      do k=0, kt+1
            do j=0, jt+1
                  thetavo(j,k) = thetao(j,k)*(1.+(0.61*qvo(j,k)))
                  thetav(j,k) = thetao(j,k)*(1.+(0.61*qvo(j,k)-qc(j,k)))
                  pio(j,k) = 100000. - (((g/(Cp*thetavo(j,k)))*(k*dk))*100.*2000.)
                  thetal(j,k) = theta(j,k) - ((La/(Cp*pio(j,k)))*qc(j,k))
            end do
      end do


      do k=1,kt
            do j=1,jt
            T = theta(j,k) * (( pio(j,k) / 100000. ) ** ( rgas / cp ))

            !Wexler's Formula
            esat = es(T)
            qvs = 0.622*(esat/(pio(j,k)-esat))
!            qvo(j,k) = qvo(j,k) + .005!RELHUM*qvs
            qv(j,k) = qvo(j,k)
            qw(j,k) = qc(j,k)+qv(j,k)
            end do
      end do

      thetao(floor(jt/2.),0) = 288.+20. !solar panel
      theta(floor(jt/2.),0) = 288.+20.

!!!!!!!!!!!!!!!!!!!!! End of initialization !!!!!!!!!!!!!!!!!!!!!!
    
      CALL BOUND

      if (animate) then
!             qc(:,kt+1) = .0007
!             qc(:,0) = 0.0

            open(71, file='av.dat', action='write',position='rewind')
            open(72, file='aw.dat', action='write',position='rewind')
            open(73, file='atheta.dat', action='write',position='rewind')
            open(74, file='api.dat', action='write',position='rewind')
            open(75, file='aqc.dat', action='write',position='rewind')
                        do k=0, kt+1
                              write(71,*) v(:,k)
                              write(72,*) w(:,k)
                              write(73,*) thetav(:,k)-thetavo(:,k)
                              write(74,*) pi(:,k)
                              write(75,*) qc(:,k)
                        end do
                  close(71)
                  close(72)
                  close(73)
                  close(74)
                  close(75)
      end if



      ITT = 1 ! itt is time step index

!     USE FORWARD SCHEME TO do first step
      A = 1.
      B = 0.
      N1 = MOD ( ITT    , 2 ) + 1
      N2 = MOD ( ITT - 1, 2 ) + 1

      CALL STEP ( N1, N2, A, B ) ! do first time step

!!!!!!! CALL ADJUST !!!!!!!!!!!

      do j = 1, jt
      do k = 1, kt
      thetaadj = theta(j,k)
      qvadj = qw(j,k)
      qcadj = 0.0
      piadj = pio(j,k)+pi(j,k)

      if (cloudtxt) then
            write(*,*) 'before: theta = ',theta(j,k), 'qw = ', qw(j,k),'qv = ',qvadj,'qc = ',qcadj!,'qvs = ',qvs
      end if
      
      CALL ADJUST (thetaadj, qvadj, qcadj, piadj)!, qvs)

      qv(j,k) = qvadj
      qc(j,k) = qcadj
      qw(j,k) = qvadj + qcadj
      theta(j,k) = thetaadj
      thetal(j,k) = thetaadj - ((La/(Cp*piadj))*qcadj)
      thetav(j,k) = theta(j,k) + (thetao(j,k)*(0.61*qv(j,k)-qc(j,k)))

      
      if (cloudtxt) then
            write(*,*) 'after:  theta = ',theta(j,k), 'qw = ', qw(j,k),'qv = ',qvadj,'qc = ',qcadj!, 'qvs = ',qvs
      end if
      end do
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      open(81, file='tv.dat') ! Open the file 'tv.dat'
      open(82, file='tw.dat') ! Open the file 'tw.dat'
      open(83, file='ttheta.dat') ! Open the file 'ttheta.dat'
      open(84, file='tpi.dat') ! Open the file 'tpi.dat'



!     ADAMS - BASHFORTH TWO - LEVEL SCHEME

      A =   3. / 2. 
      B = - 1. / 2. 

      ITTNOW = ITT  + 1

      DO ITT = ITTNOW, ITTMAX

      N1 = MOD ( ITT    , 2 ) + 1   !N1 and N2 alternnate values of 1 and 2 here
      N2 = MOD ( ITT - 1, 2 ) + 1

      CALL STEP ( N1, N2, A, B ) ! do subsequent time steps

!!!!!!! CALL ADJUST !!!!!!!!!!!

      do j = 1, jt
      do k = 1, kt
      thetaadj = theta(j,k)
      qvadj = qw(j,k)
      qcadj = 0.0
      piadj = pio(j,k)+pi(j,k)

      if (cloudtxt) then
            write(*,*) 'before: thetal = ',thetal(j,k), 'qw = ', qw(j,k),'qv = ',qvadj,'qc = ',qcadj!,'qvs = ',qvs
      end if
      
      CALL ADJUST (thetaadj, qvadj, qcadj, piadj)!, qvs)

      qv(j,k) = qvadj
      qc(j,k) = qcadj
      qw(j,k) = qvadj + qcadj
      theta(j,k) = thetaadj
      thetal(j,k) = thetaadj - ((La/(Cp*piadj))*qcadj)
      thetav(j,k) = theta(j,k) + (thetao(j,k)*(0.61*qv(j,k)-qc(j,k)))

      
      if (cloudtxt) then
            write(*,*) 'after:  thetal = ',thetal(j,k), 'qw = ', qw(j,k),'qv = ',qvadj,'qc = ',qcadj!, 'qvs = ',qvs
      end if
      end do
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! IF (mod(1.*itt-1.,1000.) .eq. 0) then
!       write(*,*) qc
! end if

      end do
	
!     END-OF-RUN OUTPUT ROUTINES GO HERE  

!            qc(:,kt+1) = .0000001
!            qc(:,0) = 0.0

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
                        write(93,*) thetav(:,k)-thetavo(:,k)
                  end do
            close(93)


           open(94, file='pi.dat') ! Open the file 'pi.dat'
                  do k=0, kt+1
                        write(94,*) pi(:,k)
                  end do
           close(94)

           open(95, file='qc.dat') ! Open file 'qc.dat'
                  do k=0, kt+1
                        write(95,*) qc(:,k)
                  end do
            close(95)

            open(96, file='thetao.dat')
                  do k=0, kt+1
                        write(96,*) thetao(:,k)
                  end do
            close(96)

            open(97, file='params.dat', action='write',position='rewind') !file for model parameters to be used in Matlab
                  write(97,*) H
                  write(97,*) dk
                  write(97,*) real(kt)
                  write(97,*) L
                  write(97,*) real(jt)
            close(97)

            open(98, file='pio.dat')
                  do k=0, kt+1
                        write(98,*) pio(:,k)
                  end do
            close(98)

           close(81)
           close(82)
           close(83)
           close(84)

           if (debug) then
                  write(*,*) '----------REALS-------'
                  write(*,*) 'dk=',dk
                  write(*,*) 'ekth=', ekth 
                  write(*,*) 'ekp=',ekp
                  write(*,*) 'ekv=', ekv
                  write(*,*) 'th0=', th0
                  write(*,*) 'Cs=', Cs
                  write(*,*) 'H=', H
                  write(*,*) 'L=', L
                  write(*,*) 'dj=', dj
                  write(*,*) 'Cp=', Cp
                  write(*,*) 'g=', g
                  write(*,*) 'tmax=', tmax
                  write(*,*) 'dt=', dt
                  write(*,*) 'delth=', delth
                  write(*,*) 'a=', a
                  write(*,*) 'b=', b
                  write(*,*) '-------------INTEGERS-----------'
                  write(*,*) 'jt=', jt
                  write(*,*) 'kt=', kt
                  write(*,*) 'jv=', jv
                  write(*,*) 'kv=', kv
                  write(*,*) 'jw=', jw
                  write(*,*) 'kw=', kw
                  write(*,*) 'jth=', jth
                  write(*,*) 'kth=', kth
                  write(*,*) 'jp=', jp
                  write(*,*) 'kp=', kp
                  write(*,*) 'ittmax=', ittmax
                  write(*,*) 'Nout=', Nout
                  write(*,*) 'j=', j
                  write(*,*) 'k=', k
                  write(*,*) 'ITTNOW=', ITTNOW
                  write(*,*) 'itt=', itt
                  write(*,*) 'N1=', N1
                  write(*,*) 'N2=', N2
           end if                                 

      write(*,*) 'End of script, saved tv, tw, ttheta, tpi, v, w, theta, and pi .dat files'



contains

      SUBROUTINE STEP ( N1, N2, A, B )

      REAL, INTENT(IN) :: A, B
      INTEGER, INTENT(IN) :: N1, N2


!     This is the entire subroutine.
      CALL RCALC ( N2 ) ! calculate forcing terms from variables at current time
      CALL AB ( N1, N2, A, B ) ! update variables using a time scheme
      CALL BOUND ! apply boundary conditions to variables

      if (mod(itt*1.,60.) .eq. 0) then

            write(81,*) sum(v)


            write(82,*) sum(w)


            write(83,*) sum(theta)


            write(84,*) sum(pi)

            if (animate) then
!                  qc(:,kt+1) = .0000001
!                  qc(:,0) = 0.0

                  open(71, file='av.dat', action='write',position='append')
                  open(72, file='aw.dat', action='write',position='append')
                  open(73, file='atheta.dat', action='write',position='append')
                  open(74, file='api.dat', action='write',position='append')
                  open(75, file='aqc.dat', action='write',position='append')
            
                        do k=0, kt+1
                              write(71,*) v(:,k)
                              write(72,*) w(:,k)
                              write(73,*) thetav(:,k)-thetavo(:,k)
                              write(74,*) pi(:,k)
                              write(75,*) qc(:,k)
                        end do
                  close(71)
                  close(72)
                  close(73)
                  close(74)
                  close(75)

                  !Stop animation output if model blows up
                  if (isnan(v(5,5))) then !checks for NaNs
                        animate = .false.
                  end if
            end if

      end if

      if (debug) then
            write(*,*) (float(itt)*100./float(ittmax)),'% Complete.'
      end if

      END SUBROUTINE STEP

      SUBROUTINE RCALC (  N2 )

      INTEGER, INTENT(IN) :: N2

!     CALCULATES FORCING TERMS FOR V(J,K), ETC.; STORES THEM  IN FV(J,K,N2), ETC.

      DO K = 1, KT
      DO J = 1, JT
      fv(J,K,N2)  = - (v(j,k)*((v(j+1,k)-v(j-1,k))/(2.*dj)))                    &
                    - (((((w(j,k)+w(j+1,k))/2.)*((v(j,k+1)-v(j,k))/dk))         &
                     +(((w(j,k-1)+w(j+1,k-1))/2.)*((v(j,k)-v(j,k-1))/dk)))/2.)  &
                    - (Cp*thetavo(j,k)*((pi(j+1,k) - pi(j,k))/dj))              &
                    + (ekv*(v(j,k+1) - (2.*v(j,k)) + v(j,k-1)) / (dk**2.))      &
                    + (ekv*(v(j+1,k) - (2.*v(j,k)) + v(j-1,k)) / (dj**2.))
      END DO
      END DO


      DO K = 1, KT
      DO J = 1, JT
      fw(j,k,N2) = - (((((v(j,k+1)+v(j,k))/2.)*((w(j+1,k)-w(j,k))/dj))              &
                        +(((v(j-1,k+1)+v(j-1,k))/2.)*((w(j,k)-w(j-1,k))/dj)))/2.)   &
            - (w(j,k)*((w(j,k+1)-w(j,k-1))/(2.*dk)))                                &
            - (Cp*(thetavo(j,k+1)+thetao(j,k))*0.5*(pi(j,k+1)-pi(j,k))/dk)          &
            + ((g/thetav(j,k))*(thetav(j,k)-thetavo(j,k)))                          &
            + (ekp*(w(j,k+1) - (2.*w(j,k)) + w(j,k-1)) / (dk**2.))                  &
            + (ekp*(w(j+1,k) - (2.*w(j,k)) + w(j-1,k)) / (dj**2.))
      END DO
      END DO


      DO K = 1, KT
      DO J = 1, JT
      ftheta(J,K,N2) =  - (((v(j,k)*((theta(j+1,k)-theta(j,k))/dj))                          &
                           +(v(j-1,k)*((theta(j,k)-theta(j-1,k))/dj)))/2.)                   &
                        + (ekth*(theta(J,k+1) - (2.*theta(J,k)) + theta(J,k-1)) / (dk**2.))   &
                        + (ekth*(theta(j+1,k) - (2.*theta(j,k)) + theta(j-1,k)) / (dj**2.))   &
                        - (0.5*(w(j,k)  *(theta(j,k+1) - theta(j,k  ))                     &
                            + (w(j,k-1)*(theta(j,k  ) - theta(j,k-1)))) / dk)
      END DO
      END DO


      DO K = 1, KT
      DO J = 1, JT
      fpi(j,k,N2) = -((Cs**2.)/(Cp*(thetavo(j,k)**2.)))                              &
                  * ((((thetavo(j,k)*v(j,k)) - (thetavo(j,k)*v(j-1,k))) / dj)      &
                  + (((thetavo(j,k+1) + thetavo(j,k))*w(j,k)) - ((thetavo(j,k) + thetavo(j,k-1))*w(j,k-1)))/(2.*dk))

      END DO
      END DO

      do k=1, kt
      do j=1,jt
      fthetal(j,k,N2) = - (((v(j,k)*((thetal(j+1,k)-thetal(j,k))/dj))                          &
                           +(v(j-1,k)*((thetal(j,k)-thetal(j-1,k))/dj)))/2.)                   &
                        + (ekth*(thetal(J,k+1) - (2.*thetal(J,k)) + thetal(J,k-1)) / (dk**2.))   &
                        + (ekth*(thetal(j+1,k) - (2.*thetal(j,k)) + thetal(j-1,k)) / (dj**2.))   &
                        - (0.5*(w(j,k)  *(thetal(j,k+1) - thetal(j,k  ))                     &
                            + (w(j,k-1)*(thetal(j,k  ) - thetal(j,k-1)))) / dk)
      end do
      end do 

      do k=1, kt
      do j=1,jt
      fqw(j,k,N2) = - (((v(j,k)*((qw(j+1,k)-qw(j,k))/dj))                          &
                           +(v(j-1,k)*((qw(j,k)-qw(j-1,k))/dj)))/2.)                   &
                        + (ekth*(qw(J,k+1) - (2.*qw(J,k)) + qw(J,k-1)) / (dk**2.))   &
                        + (ekth*(qw(j+1,k) - (2.*qw(j,k)) + qw(j-1,k)) / (dj**2.))   &
                        - (0.5*(w(j,k)  *(qw(j,k+1) - qw(j,k  ))                     &
                            + (w(j,k-1)*(qw(j,k  ) - qw(j,k-1)))) / dk)
      end do
      end do            


!     ETC (forcing for w, theta, and pi)

      END SUBROUTINE RCALC

      SUBROUTINE AB ( N1, N2, A, B )

      REAL, INTENT(IN) :: A, B
      INTEGER, INTENT(IN) :: N1, N2

!     THE FOLLOWING LOOP UPDATES V USING EITHER THE FORWARD OR THE ADAMS-BASHFORTH 
!     SCHEME DEPENDING ON THE VALUES OF A, B.  
!     SUBSCRIPT N2 OF FV ALWAYS REFERS TO THE MOST RECENTLY CALCULATED VALUES FOR FV.

      DO K = 1, KT
      DO J = 1, JT
      V(J,K)      = V(J,K)      + (DT  * (A * FV(J,K,N2)      + B * FV(J,K,N1)))
      w(j,k)      = w(j,k)      + (DT  * (A * fw(j,k,N2)      + B * fw(j,k,N1)))
      theta(J,K)  = theta(J,K)  + (DT  * (A * ftheta(J,K,N2)  + B * ftheta(J,K,N1)))
      thetal(J,K) = thetal(J,K) + (DT  * (A * fthetal(J,K,N2) + B * fthetal(J,K,N1)))
      pi(j,k)     = pi(j,k)     + (DT  * (A * fpi(j,k,N2)     + B * fpi(j,k,N1)))
      qw(j,k)     = qw(j,k)     + (DT  * (A * fqw(j,k,N2)     + B * fqw(j,k,N1)))

!            qv(j,k) = qv(j,k) + .00001!!!!!!!!!!!!!!!!!!!!!!!!!!
!            qc(j,k) = qc(j,k) + 0.007*(exp((k*dk)/(kt*dk)))

! CALL ADJUST

!       T = thetal(j,k) * (( pio(j,k) / 100000. ) ** ( rgas / cp ))

!       !Wexler's Formula
!       es = exp(   (-0.29912729e4  *((T)**(-2)))  + &
!       (-0.60170128e4  *((T)**(-1)))  + &
!             ( 0.1887643854e2*((T)**( 0)))  + & 
!             (-0.28354721e-1 *((T)**( 1)))  + &
!             ( 0.17838301e-4 *((T)**( 2)))  + &
!             (-0.84150417e-9 *((T)**( 3)))  + &
!             ( 0.44412543e-12*((T)**( 4)))  + &
!             ( 0.2858487e1*log( T)))

!       qvs = 0.622*(es/(pio(j,k)-es))



!       if (qw(j,k) .gt. qvs) then
!             qv(j,k) = qvs
!             qc(j,k) = qw(j,k) - qv(j,k)
!             theta(j,k) = thetal(j,k) + ((L/(Cp*pio(j,k))*qc(j,k)))
!       else
!             qc(j,k) = 0
!             qv(j,k) = qw(j,k)
!             theta(j,k) = thetal(j,k)
!       end if

      thetav(j,k) = theta(j,k) + thetao(j,k)*((0.61*qv(j,k))-qc(j,k))
      
                  if ((isnan(v(j,k))) .or. (isnan(w(j,k)))) then !checks for NaNs
                        write(*,*) "Model blew up at time = ",itt*dt," seconds"
                        stop
                  end if

      END DO
      END DO


!     ETC (update w, theta, and pi)

      END SUBROUTINE AB

      SUBROUTINE BOUND


!     apply boundary conditions for v

      do j = 0, jt+1 
      v(j,0) = v(j,1) ! free slip b.c.
      v(j,kt+1) = v(j,kt) ! free slip b.c.
      w(j,0) = 0. ! no vert. motion at surface
      w(j,kt) = 0.
      w(j,kt+1) = 0.
      theta(j,0) = thetao(j,0)+thetao(j,1)-theta(j,1)
      theta(j,kt+1) = thetao(j,kt+1)+thetao(j,kt)-theta(j,kt)
      thetav(j,0) = thetavo(j,0)+thetavo(j,1)-thetav(j,1)
      thetav(j,kt+1) = thetavo(j,kt+1)+thetavo(j,kt)-thetav(j,kt)
      pi(j,0) = pi(j,1)
      pi(j,kt+1) = pi(j,kt)
      qv(j,0) = qv(j,1)
      qv(j,kt+1) = qv(j,kt)
      qvo(j,0) = qvo(j,1)
      qvo(j,kt+1) = qvo(j,kt)
      qw(j,0) = qw(j,1)
      qw(j,kt+1) = qw(j,kt)
      qc(j,0) = qc(j,1)
      qc(j,kt+1) = qc(j,kt)
      thetal(j,0) = theta(j,0) - ((La/(Cp*pio(j,0)))*qc(j,0))
      thetal(j,kt+1) = theta(j,kt+1) - ((La/(Cp*pio(j,kt+1)))*qc(j,kt+1))
      end do

      do k = 1, kt
      v(0,k) = v(jt,k) ! periodic b.c.
      v(jt+1,k) = v(1,k) ! periodic b.c.
      theta(0,k) = theta(jt,k)
      theta(jt+1,k) = theta(1,k)
      pi(0,k) = pi(jt,k)
      pi(jt+1,k) = pi(1,k)
      thetav(0,k) = thetav(jt,k)
      thetav(jt+1,k) = thetav(1,k)
      thetal(0,k) = thetal(jt,k)
      thetal(jt+1,k) = thetal(1,k)
      qv(0,k) = qv(jt,k)
      qv(jt+1,k) = qv(1,k)
      qvo(0,k) = qvo(jt,k)
      qvo(jt+1,k) = qvo(1,k)
      qw(0,k) = qw(jt,k)
      qw(jt+1,k) = qw(1,k)
      qc(0,k) = qc(jt,k)
      qc(jt+1,k) = qc(1,k)
      end do

      do k=1, kt-1
      w(0,k) = w(jt,k)
      w(jt+1,k) = w(1,k)
      end do


!     ETC (apply b.c. for w, theta, and pi)

      END SUBROUTINE BOUND


      REAL FUNCTION es(T)


            REAL, INTENT(IN) :: T

      ! Wexler's formula for es(T)
            es = exp(   (-0.29912729e4  *((T)**(-2)))  + &
                  (-0.60170128e4  *((T)**(-1)))  + &
                        ( 0.1887643854e2*((T)**( 0)))  + & 
                        (-0.28354721e-1 *((T)**( 1)))  + &
                        ( 0.17838301e-4 *((T)**( 2)))  + &
                        (-0.84150417e-9 *((T)**( 3)))  + &
                        ( 0.44412543e-12*((T)**( 4)))  + &
                        ( 0.2858487e1*log( T)))

      END FUNCTION

      REAL FUNCTION DESDT(T)


            REAL TC
            REAL, INTENT(IN) :: T

            TC = T-273.15

            if (TC .lt. -50.) then
                  TC = -50.
            end if

      !     LOWE'S FORMULA FOR THE DERIVATIVE OF
      !     SATURATION VAPOR PRESSURE WITH RESPECT TO TEMPERATURE.
      !     ES IS IN PASCALS. T IS IN DEGREES KELVIN.

            DESDT = 100. * ((((( ( (-7.090245E-13*TC) + &
                                               3.532422E-10*TC)   + &
                                               1.036561E-07*TC)   + &
                                               1.215215E-05*TC)   + &
                                                 7.938054E-04*TC)       + &
                                                 2.857003E-02*TC)       + &
                                                 4.438100E-01)

      END FUNCTION
      
end program qcom
