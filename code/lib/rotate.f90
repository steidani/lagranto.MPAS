! ********************************************************************************
! * Transformation routines from gm2em                                           *
! ********************************************************************************

REAL FUNCTION LMSTOLM (PHIS, LAMS, POLPHI, POLLAM)
  !
  !**** LMSTOLM  -   FC:BERECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE FUER
  !****                 EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
  !****                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
  !****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
  !**   AUFRUF   :   LAM = LMSTOLM (PHIS, LAMS, POLPHI, POLLAM)
  !**   ENTRIES  :   KEINE
  !**   ZWECK    :   BERECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE FUER
  !**                EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
  !**                IM ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
  !**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
  !**   VERSIONS-
  !**   DATUM    :   03.05.90
  !**
  !**   EXTERNALS:   KEINE
  !**   EINGABE-
  !**   PARAMETER:   PHIS     REAL   GEOGR. BREITE DES PUNKTES IM ROT.SYS.
  !**                LAMS     REAL   GEOGR. LAENGE DES PUNKTES IM ROT.SYS.
  !**                POLPHI   REAL   WAHRE GEOGR. BREITE DES NORDPOLS
  !**                POLLAM   REAL   WAHRE GEOGR. LAENGE DES NORDPOLS
  !**   AUSGABE-
  !**   PARAMETER:   WAHRE GEOGRAPHISCHE LAENGE ALS WERT DER FUNKTION
  !**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
  !**
  !**   COMMON-
  !**   BLOECKE  :   KEINE
  !**
  !**   FEHLERBE-
  !**   HANDLUNG :   KEINE
  !**   VERFASSER:   D.MAJEWSKI

  REAL :: LAMS,PHIS,POLPHI,POLLAM

  DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /

  ZSINPOL = SIN(ZPIR18*POLPHI)
  ZCOSPOL = COS(ZPIR18*POLPHI)
  ZLAMPOL = ZPIR18*POLLAM
  ZPHIS   = ZPIR18*PHIS
  ZLAMS   = LAMS
  IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
  ZLAMS   = ZPIR18*ZLAMS

  ZARG1   = SIN(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  + &
       ZCOSPOL*           SIN(ZPHIS)) - &
       COS(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
  ZARG2   = COS(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  + &
       ZCOSPOL*           SIN(ZPHIS)) + &
       SIN(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
  IF (ABS(ZARG2).LT.1.E-30) THEN
    IF (ABS(ZARG1).LT.1.E-30) THEN
      LMSTOLM =   0.0
    ELSEIF (ZARG1.GT.0.) THEN
          LMSTOLAM =  90.0
        ELSE
          LMSTOLAM = -90.0
        ENDIF
  ELSE
    LMSTOLM = ZRPI18*ATAN2(ZARG1,ZARG2)
  ENDIF

  RETURN
END FUNCTION LMSTOLM


REAL FUNCTION PHSTOPH (PHIS, LAMS, POLPHI, POLLAM)
  !
  !**** PHSTOPH  -   FC:BERECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE FUER
  !****                 EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
  !****                 ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
  !****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
  !**   AUFRUF   :   PHI = PHSTOPH (PHIS, LAMS, POLPHI, POLLAM)
  !**   ENTRIES  :   KEINE
  !**   ZWECK    :   BERECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE FUER
  !**                EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
  !**                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
  !**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
  !**   VERSIONS-
  !**   DATUM    :   03.05.90
  !**
  !**   EXTERNALS:   KEINE
  !**   EINGABE-
  !**   PARAMETER:   PHIS     REAL   GEOGR. BREITE DES PUNKTES IM ROT.SYS.
  !**                LAMS     REAL   GEOGR. LAENGE DES PUNKTES IM ROT.SYS.
  !**                POLPHI   REAL   WAHRE GEOGR. BREITE DES NORDPOLS
  !**                POLLAM   REAL   WAHRE GEOGR. LAENGE DES NORDPOLS
  !**   AUSGABE-
  !**   PARAMETER:   WAHRE GEOGRAPHISCHE BREITE ALS WERT DER FUNKTION
  !**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
  !**
  !**   COMMON-
  !**   BLOECKE  :   KEINE
  !**
  !**   FEHLERBE-
  !**   HANDLUNG :   KEINE
  !**   VERFASSER:   D.MAJEWSKI

  REAL :: LAMS,PHIS,POLPHI,POLLAM

  DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /

  SINPOL = SIN(ZPIR18*POLPHI)
  COSPOL = COS(ZPIR18*POLPHI)
  ZPHIS  = ZPIR18*PHIS
  ZLAMS  = LAMS
  IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
  ZLAMS  = ZPIR18*ZLAMS
  ARG     = COSPOL*COS(ZPHIS)*COS(ZLAMS) + SINPOL*SIN(ZPHIS)

  PHSTOPH = ZRPI18*ASIN(ARG)

  RETURN
END FUNCTION PHSTOPH


REAL FUNCTION LMTOLMS (PHI, LAM, POLPHI, POLLAM)
  !
  !%Z% Modul %M%, V%I% vom %G%, extrahiert am %H%
  !
  !**** LMTOLMS  -   FC:UMRECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE LAM
  !****                 AUF EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
  !****                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
  !****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
  !**   AUFRUF   :   LAM = LMTOLMS (PHI, LAM, POLPHI, POLLAM)
  !**   ENTRIES  :   KEINE
  !**   ZWECK    :   UMRECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE LAM AUF
  !**                EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
  !**                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
  !**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
  !**   VERSIONS-
  !**   DATUM    :   03.05.90
  !**
  !**   EXTERNALS:   KEINE
  !**   EINGABE-
  !**   PARAMETER:   PHI    REAL BREITE DES PUNKTES IM GEOGR. SYSTEM
  !**                LAM    REAL LAENGE DES PUNKTES IM GEOGR. SYSTEM
  !**                POLPHI REAL GEOGR.BREITE DES N-POLS DES ROT. SYSTEMS
  !**                POLLAM REAL GEOGR.LAENGE DES N-POLS DES ROT. SYSTEMS
  !**   AUSGABE-
  !**   PARAMETER:   WAHRE GEOGRAPHISCHE LAENGE ALS WERT DER FUNKTION
  !**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
  !**
  !**   COMMON-
  !**   BLOECKE  :   KEINE
  !**
  !**   FEHLERBE-
  !**   HANDLUNG :   KEINE
  !**   VERFASSER:   G. DE MORSIER

  REAL :: LAM,PHI,POLPHI,POLLAM

  DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /

  ZSINPOL = SIN(ZPIR18*POLPHI)
  ZCOSPOL = COS(ZPIR18*POLPHI)
  ZLAMPOL =     ZPIR18*POLLAM
  ZPHI    =     ZPIR18*PHI
  ZLAM    = LAM
  IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
  ZLAM    = ZPIR18*ZLAM

  ZARG1   = - SIN(ZLAM-ZLAMPOL)*COS(ZPHI)
  ZARG2   = - ZSINPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL)+ZCOSPOL*SIN(ZPHI)
  IF (ABS(ZARG2).LT.1.E-30) THEN
    IF (ABS(ZARG1).LT.1.E-30) THEN
      LMTOLMS =   0.0
    ELSEIF (ZARG1.GT.0.) THEN
          LMTOLMS =  90.0
        ELSE
          LMTOLMS = -90.0
        ENDIF
  ELSE
    LMTOLMS = ZRPI18*ATAN2(ZARG1,ZARG2)
  ENDIF

  RETURN
END FUNCTION LMTOLMS


REAL FUNCTION PHTOPHS (PHI, LAM, POLPHI, POLLAM)
  !
  !%Z% Modul %M%, V%I% vom %G%, extrahiert am %H%
  !
  !**** PHTOPHS  -   FC:UMRECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE PHI
  !****                 AUF EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
  !****                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
  !****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
  !**   AUFRUF   :   PHI = PHTOPHS (PHI, LAM, POLPHI, POLLAM)
  !**   ENTRIES  :   KEINE
  !**   ZWECK    :   UMRECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE PHI AUF
  !**                EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
  !**                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
  !**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
  !**   VERSIONS-
  !**   DATUM    :   03.05.90
  !**
  !**   EXTERNALS:   KEINE
  !**   EINGABE-
  !**   PARAMETER:   PHI    REAL BREITE DES PUNKTES IM GEOGR. SYSTEM
  !**                LAM    REAL LAENGE DES PUNKTES IM GEOGR. SYSTEM
  !**                POLPHI REAL GEOGR.BREITE DES N-POLS DES ROT. SYSTEMS
  !**                POLLAM REAL GEOGR.LAENGE DES N-POLS DES ROT. SYSTEMS
  !**   AUSGABE-
  !**   PARAMETER:   ROTIERTE BREITE PHIS ALS WERT DER FUNKTION
  !**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
  !**
  !**   COMMON-
  !**   BLOECKE  :   KEINE
  !**
  !**   FEHLERBE-
  !**   HANDLUNG :   KEINE
  !**   VERFASSER:   G. DE MORSIER

  REAL :: LAM,PHI,POLPHI,POLLAM

  DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /

  ZSINPOL = SIN(ZPIR18*POLPHI)
  ZCOSPOL = COS(ZPIR18*POLPHI)
  ZLAMPOL = ZPIR18*POLLAM
  ZPHI    = ZPIR18*PHI
  ZLAM    = LAM
  IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
  ZLAM    = ZPIR18*ZLAM
  ZARG    = ZCOSPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL) + ZSINPOL*SIN(ZPHI)

  PHTOPHS = ZRPI18*ASIN(ZARG)

  RETURN
END FUNCTION PHTOPHS


SUBROUTINE uv2uvrot(u, v, rlat, rlon, pollat, pollon, urot, vrot)

  ! Description:
  ! This routine converts the wind components u and v from the real
  ! geographical system to the rotated system.
  !
  ! Method:
  ! Transformation formulas for converting between these two systems.

  real :: u, v
  real :: rlat, rlon
  real :: pollat, pollon
  real :: urot, vrot
  real :: zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm
  real      :: zrpi18 = 57.2957795
  real      :: zpir18 = 0.0174532925

  zsinpol = SIN(pollat * zpir18)
  zcospol = COS(pollat * zpir18)
  zlonp   = (pollon-rlon) * zpir18
  zlat    =         rlat  * zpir18
  zarg1   = zcospol*SIN(zlonp)
  zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
  znorm   = 1./SQRT( zarg1**2 + zarg2**2 )
  urot   =  u*zarg2*znorm - v*zarg1*znorm
  vrot   =  u*zarg1*znorm + v*zarg2*znorm

END SUBROUTINE uv2uvrot


SUBROUTINE uvrot2uv (urot, vrot, rlat, rlon, pollat, pollon, u, v)

  ! Description:
  ! This routine converts the wind components u and v from the rotated system
  ! to the real geographical system.
  !
  ! Method:
  ! Transformation formulas for converting between these two systems.

  integer :: n
  real :: u, v
  real :: rlat, rlon
  real :: pollat, pollon
  real :: urot, vrot
  real :: zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm
  integer :: i
  real      :: zrpi18 = 57.2957795
  real      :: zpir18 = 0.0174532925

  zsinpol = SIN(pollat * zpir18)
  zcospol = COS(pollat * zpir18)
  zlonp   = (pollon-rlon) * zpir18
  zlat    =         rlat  * zpir18
  zarg1   = zcospol*SIN(zlonp)
  zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
  znorm   = 1./SQRT(zarg1**2 + zarg2**2)
  u       =   urot*zarg2*znorm + vrot*zarg1*znorm
  v       = - urot*zarg1*znorm + vrot*zarg2*znorm

END SUBROUTINE uvrot2uv

! ********************************************************************************
! * Transformation routines from <utilities.f90>                                 *
! ********************************************************************************

FUNCTION  phirot2phi ( phirot, rlarot, polphi, pollam, polgam )

!------------------------------------------------------------------------------
!
! Description:
!   This function converts phi from one rotated system to phi in another
!   system. If the optional argument polgam is present, the other system
!   can also be a rotated one, where polgam is the angle between the two
!   north poles.
!   If polgam is not present, the other system is the real geographical
!   system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

! Parameter list:
REAL, INTENT (IN)      ::        &
  polphi,   & ! latitude of the rotated north pole
  pollam,   & ! longitude of the rotated north pole
  phirot,   & ! latitude in the rotated system
  rlarot      ! longitude in the rotated system

REAL, INTENT (IN)      ::        &
  polgam      ! angle between the north poles of the systems

REAL                   ::        &
  phirot2phi  ! latitude in the geographical system

! Local variables
REAL                   ::        &
  zsinpol, zcospol, zphis, zrlas, zarg, zgam

REAL, PARAMETER        ::        &
  zrpi18 = 57.2957795,                  &
  zpir18 = 0.0174532925

!------------------------------------------------------------------------------

! Begin function phirot2phi

  zsinpol     = SIN (zpir18 * polphi)
  zcospol     = COS (zpir18 * polphi)

  zphis       = zpir18 * phirot
  IF (rlarot > 180.0) THEN
    zrlas = rlarot - 360.0
  ELSE
    zrlas = rlarot
  ENDIF
  zrlas       = zpir18 * zrlas

  IF (polgam /= 0.0) THEN
    zgam  = zpir18 * polgam
    zarg  = zsinpol*SIN (zphis) +                                           &
        zcospol*COS(zphis) * ( COS(zrlas)*COS(zgam) - SIN(zgam)*SIN(zrlas) )
  ELSE
    zarg  = zcospol * COS (zphis) * COS (zrlas) + zsinpol * SIN (zphis)
  ENDIF

  phirot2phi  = zrpi18 * ASIN (zarg)

END FUNCTION phirot2phi


FUNCTION  phi2phirot ( phi, rla, polphi, pollam )

!------------------------------------------------------------------------------
! Description:
!   This routine converts phi from the real geographical system to phi
!   in the rotated system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------
! Parameter list:
REAL, INTENT (IN)      ::        &
  polphi,  & ! latitude of the rotated north pole
  pollam,  & ! longitude of the rotated north pole
  phi,     & ! latitude in the geographical system
  rla        ! longitude in the geographical system

REAL                   ::        &
  phi2phirot ! longitude in the rotated system

! Local variables
REAL                       ::    &
  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

REAL, PARAMETER            ::    &
  zrpi18 = 57.2957795,                  & !
  zpir18 = 0.0174532925

!------------------------------------------------------------------------------

! Begin function phi2phirot

  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0) THEN
    zrla1  = rla - 360.0
  ELSE
    zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1

  zarg1    = SIN (zphi) * zsinpol
  zarg2    = COS (zphi) * zcospol * COS (zrla - zlampol)

  phi2phirot = zrpi18 * ASIN (zarg1 + zarg2)

END FUNCTION phi2phirot

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

FUNCTION  rlarot2rla (phirot, rlarot, polphi, pollam, polgam)

!------------------------------------------------------------------------------
!
! Description:
!   This function converts lambda from one rotated system to lambda in another
!   system. If the optional argument polgam is present, the other system
!   can also be a rotated one, where polgam is the angle between the two
!   north poles.
!   If polgam is not present, the other system is the real geographical
!   system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
! Modules used:    NONE
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

! Parameter list:
REAL, INTENT (IN)      ::        &
  polphi,   & ! latitude of the rotated north pole
  pollam,   & ! longitude of the rotated north pole
  phirot,   & ! latitude in the rotated system
  rlarot      ! longitude in the rotated system

REAL, INTENT (IN)      ::        &
  polgam      ! angle between the north poles of the systems

REAL                   ::        &
  rlarot2rla  ! latitude in the geographical system

! Local variables
REAL                   ::        &
  zsinpol, zcospol, zlampol, zphis, zrlas, zarg1, zarg2, zgam

REAL, PARAMETER        ::        &
  zrpi18 = 57.2957795,                  & !
  zpir18 = 0.0174532925

!------------------------------------------------------------------------------

! Begin function rlarot2rla

  zsinpol = SIN (zpir18 * polphi)
  zcospol = COS (zpir18 * polphi)

  zlampol = zpir18 * pollam
  zphis   = zpir18 * phirot
  IF (rlarot > 180.0) THEN
    zrlas = rlarot - 360.0
  ELSE
    zrlas = rlarot
  ENDIF
  zrlas   = zpir18 * zrlas

  IF (polgam /= 0.0) THEN
    zgam    = zpir18 * polgam
    zarg1   = SIN (zlampol) *                                                &
      (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
       + zcospol * SIN(zphis))                                               &
    - COS (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))

    zarg2   = COS (zlampol) *                                                &
      (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
       + zcospol * SIN(zphis))                                               &
    + SIN (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))
  ELSE
    zarg1   = SIN (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
                                zcospol *              SIN(zphis)) -    &
              COS (zlampol) *             SIN(zrlas) * COS(zphis)
    zarg2   = COS (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
                                zcospol *              SIN(zphis)) +   &
              SIN (zlampol) *             SIN(zrlas) * COS(zphis)
  ENDIF

  IF (zarg2 == 0.0) zarg2 = 1.0E-20

  rlarot2rla = zrpi18 * ATAN2(zarg1,zarg2)

END FUNCTION rlarot2rla

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

FUNCTION  rla2rlarot ( phi, rla, polphi, pollam, polgam )

!------------------------------------------------------------------------------
!
! Description:
!   This routine converts lambda from the real geographical system to lambda
!   in the rotated system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------
!
! Parameter list:
REAL, INTENT (IN)      ::        &
  polphi,  & ! latitude of the rotated north pole
  pollam,  & ! longitude of the rotated north pole
  phi,     & ! latitude in geographical system
  rla        ! longitude in geographical system

REAL, INTENT (IN)      ::        &
  polgam      ! angle between the north poles of the systems

REAL                   ::        &
  rla2rlarot ! latitude in the the rotated system

! Local variables
REAL                       ::    &
  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

REAL, PARAMETER            ::    &
  zrpi18 = 57.2957795,                  & !
  zpir18 = 0.0174532925

!------------------------------------------------------------------------------

! Begin function rla2rlarot

  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0) THEN
    zrla1  = rla - 360.0
  ELSE
    zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1

  zarg1    = - SIN (zrla-zlampol) * COS(zphi)
  zarg2    = - zsinpol * COS(zphi) * COS(zrla-zlampol) + zcospol * SIN(zphi)

  IF (zarg2 == 0.0) zarg2 = 1.0E-20

  rla2rlarot = zrpi18 * ATAN2 (zarg1,zarg2)

  IF (polgam /= 0.0 ) THEN
    rla2rlarot = polgam + rla2rlarot
    IF (rla2rlarot > 180.) rla2rlarot = rla2rlarot -360.
  ENDIF

END FUNCTION rla2rlarot

