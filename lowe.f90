REAL FUNCTION DESDT (T)

implicit none

	REAL TC
	REAL, INTENT(IN) :: T
	REAL, dimension (7) :: COEF

	COEF(1:7) = (/4.438100E-01,2.857003E-02,7.938054E-04,1.215215E-05, &
					1.036561E-07,3.532422E-10,-7.090245E-13/)


	TC = T-273.15

	if (TC .lt. -50.) then
		TC = -50.
	end if

!     LOWE'S FORMULA FOR THE DERIVATIVE OF
!     SATURATION VAPOR PRESSURE WITH RESPECT TO TEMPERATURE.
!     ES IS IN PASCALS. T IS IN DEGREES KELVIN.

 	DESDT = 100. * ((((((COEF(7))   	* &
 						(TC + COEF(6))  * &
 						(TC + COEF(5))) * &
 						(TC + COEF(4))) * &
						(TC + COEF(3))) * &
						(TC + COEF(2))) * &
						(TC * COEF(1)))

END FUNCTION