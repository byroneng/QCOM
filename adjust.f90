SUBROUTINE adjust(THETAa, QVa, QCa, PBARa)

! This subroutine performs a saturation adjustment.
! Input values:
! pbar - pressure (Pa)

! In/Out values:
! theta - potential temperature (K)
! qv - water vapor mixing ratio (kg/kg)
! qc - liquid water mixing ratio (kg/kg)

implicit none

	REAL, INTENT(IN) :: PBARa
	REAL, INTENT(INOUT) :: THETAa, QVa, QCa
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

PIBAR = (PBARa/PZERO)**(RGAS/CP)
GAMMA = HLF / (CP*PIBAR)

THSTAR = THETAa
QVSTAR = QVa

TSTAR = THSTAR * PIBAR
ES1 = ES(TSTAR)

ALPHA = DESDT(TSTAR) * EPS * PIBAR * PBARa / (PBARa - ES1)**2
THFAC = GAMMA / (1. + GAMMA*ALPHA)

QVSAT = EPS / (PBARa - ES1) * ES1
THETAa = THSTAR + THFAC * (QVSTAR - QVSAT)
QW = QVSTAR + QCa !old qv

QVa = QVSAT + ALPHA * (THETAa-THSTAR) !new qv

QCa = QW - QVa !new qv

if (QCa .LT. 0.) then
	QCa = 0.
	QVa = QW
	THETAa = THSTAR + GAMMA * (QVSTAR-QVa)
end if

RETURN

END SUBROUTINE