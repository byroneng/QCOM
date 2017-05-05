REAL FUNCTION DESDT (T)

implicit none

	REAL TC
	REAL, INTENT(IN) :: T

	TC = T-273.15

	if (TC .lt. -50.) then
		TC = -50.
	end if

!     LOWE'S FORMULA FOR THE DERIVATIVE OF
!     SATURATION VAPOR PRESSURE WITH RESPECT TO TEMPERATURE.
!     ES IS IN PASCALS. T IS IN DEGREES KELVIN.

 	DESDT = 100. * ((((( ( (-7.090245E-13*TC)	+ &
 						     3.532422E-10*TC) 	+ &
 						     1.036561E-07*TC) 	+ &
 						     1.215215E-05*TC) 	+ &
							 7.938054E-04*TC) 	+ &
							 2.857003E-02*TC) 	+ &
							 4.438100E-01)

END FUNCTION