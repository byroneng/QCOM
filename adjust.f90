SUBROUTINE adjust(THETA, QV, QC, PBAR)

! This subroutine performs a saturation adjustment.
! Input values:
! pbar - pressure (Pa)

! In/Out values:
! theta - potential temperature (K)
! qv - water vapor mixing ratio (kg/kg)
! qc - liquid water mixing ratio (kg/kg)

implicit none

	REAL, INTENT(IN) :: PBAR
	REAL, INTENT(INOUT) :: THETA, QV, QC
	REAL THSTAR, QVSTAR
	REAL HLF, CP, RGAS, PZERO, PIBAR
	REAL GAMMA, TSTAR, EPS, ALPHA
	REAL QVSAT, ES1, THFAC, QW
	REAL ES, DESDT


! Constants
HLF = 2.500E+06
CP = 1004.5
RGAS = 287.04
PZERO = 100000.
EPS = 0.622

PIBAR = (PBAR/PZERO)**(RGAS/CP)
GAMMA = HLF / (CP*PIBAR)

THSTAR = THETA
QVSTAR = QV

TSTAR = THSTAR * PIBAR
ES1 = ES(TSTAR)

ALPHA = DESDT(TSTAR) * EPS * PIBAR * PBAR / (PBAR - ES1)**2
THFAC = GAMMA / (1. + GAMMA*ALPHA)

QVSAT = EPS / (PBAR - ES1) * ES1
THETA = THSTAR + THFAC * (QVSTAR - QVSAT)
QV = QVSAT + ALPHA * (THETA-THSTAR)

QW = QV + QC
QC = QW - QV

if (QC .LT. 0.) then
	QC = 0.
	QV = QW
	THETA = THSTAR + GAMMA * (QVSTAR-QV)
end if

RETURN

END SUBROUTINE