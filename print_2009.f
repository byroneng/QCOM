      SUBROUTINE PRINT ( F, JDIM, J3, NJ, K1, NK, TITLE )
      DIMENSION F(JDIM,1), IA(33), TITLE(2)
c      
c     removed carriage control characters -- 19 Feb 2009    
c
c     Prints a 2D array in scaled integer format; 33 columns at a time
c
c     Arguments:
c
c     F      2D array to print
c     JDIM   first dimension of F
c     J3     first column of F to print
c     NJ     last column of F to print
c     K1     first row of F to print (appears at bottom)
c     NK     last row of F to print (appears at top)
c     TITLE  8-character title string
C
C     REVISED 9 JULY 1982 : ADDITIONAL ARGUMENTS ( J3, K1 )
C
      J2 = J3 - 1
    1 J1 = J2 + 1
      NJP = AMIN0 ( 33, NJ - J1 + 1 )
      J2 = NJP + J1 - 1
C
      CALL PRINT1 ( F, IA, JDIM, NJP, J1, J2, K1, NK, TITLE )
C
      IF ( J2 .LT. NJ ) GO TO 1
C
      RETURN
      END
      SUBROUTINE PRINT1 ( A, IA, NJ, NJP, J1, J2, K1, K2, TITLE )
      DIMENSION A (NJ, 1), IA (NJP)
      REAL *8 TITLE
C
C     *** REVISED 13 DEC 1985 ***
C
C     *** FIND X = MAX ( ABS ( A ) )
C
    9 XMIN =   1.E12
      XMAX = - 1.E12
      DO 10 K = K1, K2
      DO 10 J = J1, J2
      XMIN = AMIN1 ( XMIN , A(J,K) )
      XMAX = AMAX1 ( XMAX , A(J,K) )
   10 CONTINUE
      X = AMAX1 ( ABS ( XMIN ), XMAX )
C
      IF(X.GT.0.) GO TO 4
C
      WRITE(6,30) TITLE
   30 FORMAT('ZERO FIELD FOR  ',A8)
      RETURN
C
C        *** FIND SCALE SUCH THAT 100 LT SCALE * X LT 1000
C        *** AND SCALE = 10 ** I, WHERE I IS AN INTEGER.
C
    4 Y = ALOG10 ( X )
      IF ( Y . LT . 0. ) GO TO 5
      IEXP = INT ( Y )
      GO TO 6
    5 IEXP = INT ( Y ) - 1
    6 SCALE = 10.** ( - IEXP + 2 )
      IEXPU = IEXP - 2
C
      IMIN = SCALE * XMIN + SIGN ( 0.5, XMIN + X * 1.E-06 )
      IMAX = SCALE * XMAX + SIGN ( 0.5, XMAX + X * 1.E-06 )
C
C        *** WRITE SCALE AND NAME OF FIELD
C
      WRITE (6,20) TITLE, IEXPU, IMIN, IMAX
   20 FORMAT (A8,5X,'( UNITS ARE 10 **',I4,' )',
     $ 5X,'EXTREMES:',2I5)
C        *** SCALE A, STORE INTO INTEGER ARRAY AND WRITE IT OUT
C
      NKP = K2 - K1 + 1
      DO 3 KK = 1, NKP
      K = K2 - KK + 1
      DO 1 JJ = 1, NJP
      J = J1 + JJ - 1
    1 IA(JJ) = SCALE * A(J,K) + SIGN ( 0.5, A(J,K) + X * 1.E-06 )
      WRITE (6,2) IA
    2 FORMAT (33I4)
    3 CONTINUE
C
      RETURN
      END
