REAL FUNCTION es (T)

implicit none

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