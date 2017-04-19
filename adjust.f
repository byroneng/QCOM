      SUBROUTINE ADJUST ( theta, qv, qc, PBAR )
c     SUBROUTINE ADJUST ( TH, QV, QC, PBAR, qvs )
c
c     Performs an isobaric moist adiabatic adjustment.
c     The final state is either subsaturated with no liquid water
c     present, or exactly saturated with liquid water present.
c
c     This version iterates to obtain the adjustment, so accurate
c     guesses for TH and QV are not needed. For use in situations
c     where TH and QV will be successively updated, as in a parcel
c     model, set ITTMAX=1.
c
c     Units: SI (MKS)
c
c     Input --
c     TH: potential temperature, theta^* (K)
c     QV: mixing ratio of water vapor, q_v^* (kg/kg)
c     QC: mixing ratio of liquid water, q_c^* (kg/kg)
c     PBAR: pressure, p (Pa)
c
c     Output --
c     TH: adjusted potential temperature, theta^{n+1} (K)
c     QV: adjusted mixing ratio of water vapor, q_v^{n+1} (kg/kg)
c     QC: adjusted mixing ratio of liquid water, q_c^{n+1} (kg/kg)
c     QVS: saturation mixing ratio based on TH and PBAR
c
c     29 Nov 91 --
c     Exact Qsat now used in ALPHA and QSAT.
c     Corrected QVS1 for QC=0.
c     5 Feb 92 --
c     PIBAR calculated within subroutine now.
c     6 Feb 92 --
c     GAMMA corrected!
c
      REAL * 8 THSTAR, TH1, QVSTAR, QV1
c      REAL * 8 QV2
      DATA HLF/2.500E+06/,CP/1004.5/
      DATA rgas/287.04/,pzero/100000./
      DATA ITTMAX/1/,DTCRIT/0.001/
C
C     PBAR IS A ( HYDROSTATIC ) REFERENCE PRESSURE FIELD.
C
      itt = 1
      pibar = ( pbar / pzero ) ** ( rgas / cp )
      GAMMA = HLF / ( CP * PIBAR )
C
      THSTAR = theta
      QVSTAR = qv
C
   30 Tstar = THSTAR * PIBAR
      es1 = es ( Tstar )
C
      ALPHA = DESDT ( Tstar ) * 0.622 * PIBAR
     $        * pbar / ( pbar - es1 ) **2
      THFAC = GAMMA / ( 1. + GAMMA * ALPHA )
C
      QVSAT = 0.622 / ( PBAR - es1 ) * ES1
      TH1 = THSTAR + THFAC * ( QVSTAR - QVSAT )
      QV1 = QVSAT + ALPHA * ( TH1 - THSTAR )
c     statement below gives same result as one above
c     QV1 = QVSTAR - ( TH1 - THSTAR ) / GAMMA
      QW1 = qv + qc
      QC1 = QW1 - QV1
C
      QVS1 = QV1
C
      IF ( QC1 .LT. 0. ) then
        QC1 = 0.
        QV1 = QW1
        TH1 = THSTAR + GAMMA * ( QVSTAR - QV1 )
        QVS1 = QVSAT + alpha * ( th1 - thstar )
      ENDIF
C
      DT = ( TH1 - THSTAR ) * PIBAR
C
      IF ( ABS( DT ) .LT. DTCRIT .OR. ITT .EQ. ITTMAX ) then
        theta = TH1
        qv = QV1
        qc = QC1
        QVS = QVS1
        return
      endif
C
      THSTAR = TH1
      QVSTAR = QV1
C
      ITT = ITT + 1
C
      GO TO 30
C
      END
      FUNCTION ES ( T )
C
C     LOWE'S FORMULA FOR SATURATION VAPOR PRESSURE ( PA ).
C     T IS IN DEGREES KELVIN.
C
      DIMENSION A(7)
      DATA A/6.107800,
     $       4.436519E-01,
     $       1.428946E-02,
     $       2.650648E-04,
     $       3.031240E-06,
     $       2.034081E-08,
     $       6.136821E-11/
C
      TC = T - 273.16
C
      IF ( TC .LT. - 50. ) TC = - 50.
C
      X = A(7)
C
      DO 1 J = 1,6
    1 X = X * TC + A(7-J)
C
      ES = X * 100.
C
      RETURN
      END
      FUNCTION DESDT ( T )
C
C     LOWE'S FORMULA FOR THE DERIVATIVE OF
C     SATURATION VAPOR PRESSURE WITH RESPECT TO TEMPERATURE.
C     ES IS IN PASCALS. T IS IN DEGREES KELVIN.
C
      DIMENSION A(7)
      DATA A/4.438100E-01,
     $       2.857003E-02,
     $       7.938054E-04,
     $       1.215215E-05,
     $       1.036561E-07,
     $       3.532422E-10,
     $     - 7.090245E-13/
C
      TC = T - 273.16
C
      IF ( TC .LT. - 50. ) TC = - 50.
C
      X = A(7)
C
      DO 1 J = 1,6
    1 X = X * TC + A(7-J)
C
      DESDT = X * 100.
C
      RETURN
      END



