! $RCSfile: utilities.f90,v $
! $Revision: 4.11 $ $Date: 2009/11/30 14:29:09 $
!+ Source module for utility routines
!==============================================================================

MODULE  utilities

!==============================================================================
!
! Description:
!   This module provides service utilities for the model. All routines are 
!   written in a manner that also other models can use it. That means:
!     - no routine uses other modules, except the declarations for the 
!       KIND-type parameter; the data access is by parameter list only
!     - no routine allocates dynamic memory; work space needed is
!       provided via the parameter list
!     - no derived data types are used
!
!   Routines (module procedures) currently contained:
!
!     - convert_month:
!       Converts a 3-character string abbreviation of a month into the number
!       of the month or vice versa.
!
!     - dfilt4:
!       Digital filter of length 4
!
!     - dfilt8:
!       Digital filter of length 8
!
!     - dolph:
!       Calculates the Dolph-Chebyshev window for the initialization
!
!     - elapsed_time:
!       Returns the elapsed wall-clock time in seconds since the last call.
!       On the first call the variables are only initialized. If no system
!       clock is present, an error-value will be returned
!
!     - get_utc_date:
!       Calculates the actual date using the date of the forecast-start and 
!       the number of timesteps performed.
!
!     - horizontal_filtering
!       horizontal filtering (at the moment especially for the pressure deviation)
!
!     - phirot2phi:
!       Converts phi from the rotated system to phi in the real
!       geographical system.
!
!     - phi2phirot:
!       Converts phi from the real geographical system to phi
!       in the rotated system.
!
!     - rlarot2rla:
!       Converts lambda from the rotated system to lambda in the real
!       geographical system.
!
!     - rla2rlarot:
!       Converts lambda from the real geographical system to lambda 
!       in the rotated system.
!
!     - sleve_split_oro
!       Decomposes a given topography field in a large-scale and a small-scale
!       part according to the definition of the SLEVE coordinate
!
!     - smoother:
!       Smoothes a 2-D field by applying digital filters
!
!     - tautsp:
!       Computes tension splines
!
!     - tautsp2D:
!       Computes tension splines for several columns
!
!     - to_upper:
!       Converts alphabetic characters from lower to upper case.
!
!     - uvrot2uv:
!       Converts the wind components u and v from the rotated system
!       to the real geographical system.
!
!     - uvrot2uv_vec:
!       the same as above, but for a whole 2D field (in vectorized form).
!
!     - uv2uvrot:
!       Converts the wind components u and v from the real geographical
!       system to the rotated system.
!
!     - uv2uvrot_vec:
!       the same as above, but for a whole 2D field (in vectorized form).
!
!     - uv2df:
!       Converts the wind components u and v to wind direction and speed.
!
!     - uv2df_vec:
!       the same as above, but for a whole 2D field (in vectorized form).
!
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.2        1998/03/30 Ulrich Schaettler
!  Introduction of subroutine dolph used during the initialization
! 1.9        1998/09/16 Guenther Doms
!  Introduction of a smoothing routine 'smoother' which uses digital
!  filters 'dfilt4' and 'dfilt8'.
! 1.10       1998/09/29 Ulrich Schaettler
!  Routine remark eliminated and put to parallel_utilities.
!  Routines uv2uvrot and uv2df introduced
! 1.16       1998/11/02 Guenther Doms
!  Correction of filter processing in routine 'smoother'.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adaptations to use this module also in GME2LM
! 1.32       1999/08/24 Guenther Doms
!  some _ireals declarations added.
! 2.8        2001/07/06 Ulrich Schaettler
!  Added new subroutines tautsp2D, uv2uvrot_vec and uvrot2uv_vec for 
!  vectorization
! 2.14       2002/02/15 Ulrich Schaettler
!  Correction and adaptations in tautsp2D (analogous to GME2LM)
!  Added new subroutine dc_topo for the SLEVE coordinate
! 2.17       2002/05/08 Ulrich Schaettler
!  Modifications for performing the filtering in irealgrib-format
! 2.18       2002/07/16 Guenther Doms
!  Corrections for the rotation of the wind components from or to the
!  geographical coordinate system.
! 3.3        2003/04/22 Christoph Schraff
!  Introduction of subroutines 'convert_month' and 'to_upper' (for GPS data).
! 3.6        2003/12/11 Ulrich Schaettler
!  Eliminated Subroutine istringlen (use F90 intrinsic LEN_TRIM instead)
! 3.13       2004/12/03 Ulrich Schaettler
!  Eliminated dependency on data_io (put irealgrib to data_parameters)
!  New SR horizontal_filtering (from INT2LM);
!  Renamed SR dc_topo to sleve_split_oro
! 3.14       2005/01/25 Ulrich Schaettler
!  New filter routine smooth9 for new type of Rayleigh damping (Lucio Torrisi)
!  Changes in horizontal_filtering (Jochen Foerstner)
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL
! 3.16       2005/07/22 Ulrich Schaettler
!  Bug correction in the call to intrinsic function REAL
! 3.18       2006/03/03 Ulrich Schaettler
!  Introduced idouble/isingle as KIND parameters instead of ireals/irealgrib
!  in the generic formulation of some routines (dfilt4, dfilt8, smoother)
!  Changed get_utc_date to include also a climatological year with 360 days
! 3.21       2006/12/04 Burkhardt Rockel, Lucio Torrisi, Jochen Foerstner
!  Added polgam in transformation function rla2rlarot
!  polgam is not used as optional parameter any more
!  Some adaptations in smooth9 for itype_spubc=2
!  Some modifications in horizontal_filtering
!  Function uv2df_vec introduced. (C. Schraff)
! V3_23        2007/03/30 Ulrich Schaettler
!  Declared some constant variables as parameters to allow inlining on
!  some platforms
!  Changed computation of acthour in get_utc_date
! V3_24        2007/04/26 Ulrich Schaettler
!  Bug correction in computation of acthour in get_utc_date
! V4_1         2007/12/04 Ulrich Schaettler
!  Introduced parameter myid to sleve_split_oro (is called from all PEs in INT2LM)
! V4_4         2008/07/16 Ulrich Schaettler
!  Adapted a debug printout in SR tautsp2D
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
!  Vectorized SR horizontal_filtering by changing some loops
!  Treatment of very small values for spline interpolation in tautsp2D
! V4_8         2009/02/16 Ulrich Schaettler (Andy Dobler)
!  Corrected leap year calculation for centuries in Gregorian calendar
! @VERSION@    @DATE@     <Your name>
!  <Modification comments>         
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
USE data_parameters , ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers, & ! KIND-type parameter for standard integer variables
    irealgrib, & ! KIND-type parameter for real variables in the grib library
    idouble,   & ! KIND-type parameter for double precision real variables
    isingle      ! KIND-type parameter for single precision real variables

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Interface Blocks
INTERFACE smoother
  MODULE PROCEDURE                        &
    smoother_double,                      &
    smoother_single
END INTERFACE

INTERFACE dfilt4
  MODULE PROCEDURE                        &
    dfilt4_double,                        &
    dfilt4_single
END INTERFACE

INTERFACE dfilt8
  MODULE PROCEDURE                        &
    dfilt8_double,                        &
    dfilt8_single
END INTERFACE

!==============================================================================

CONTAINS

!==============================================================================

SUBROUTINE convert_month ( MonthString, MonthNumber, ind )

!-------------------------------------------------------------------------------
!
! Description:
!   Convert 3-chr Month string to number (ind > 0) or vice-versa (ind <= 0).
!     If ind > 0 and input string not valid, MonthNumber will be 0.
!     If ind <=0 and input month number invalid, MonthString will be 'XXX'.
! Method:
!     Uses subroutine 'to_upper'.
!-------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  CHARACTER (LEN=3)       , INTENT(INOUT) :: MonthString  ! 3-chr Month name
  INTEGER (KIND=iintegers), INTENT(INOUT) :: MonthNumber  ! Month number (1-12)
  INTEGER (KIND=iintegers), INTENT(IN)    :: ind          ! > 0 : chr --> no.
                                                          ! <=0 : no  --> chr

! Local parameters:
! ----------------
  CHARACTER (LEN=36), PARAMETER ::  &
    MonthNames = "JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC"

! Local variables
! ---------------
  CHARACTER (LEN=3)             :: Month
  INTEGER (KIND=iintegers)      :: idx
!
!------------ End of header ----------------------------------------------------

  IF ( ind > 0 ) THEN

! ----- String to number -----

    Month = MonthString
    CALL to_upper ( Month )
    idx = INDEX ( MonthNames, Month )
    IF ( MOD ( idx-1, 3 ) == 0 ) THEN
      MonthNumber = ( idx + 2 ) / 3
    ELSE
      MonthNumber = 0
    END IF

  ELSE

! ----- Number to string -----

    IF ( MonthNumber >= 1 .AND. &
         MonthNumber <= 12 ) THEN
      idx = MonthNumber * 3 - 2
      MonthString = MonthNames(idx:idx+2)
    ELSE
      MonthString = "XXX"
    END IF

  END IF

END SUBROUTINE convert_month

!------------------------------------------------------------------------------

!==============================================================================
!==============================================================================
!+ Defines all subroutines for the generic routine dfilt4
!------------------------------------------------------------------------------
!
!  SUBROUTINE dfilt4 (fin, idim, fhelp, fout, nfilt)
!
!------------------------------------------------------------------------------
!
! Description:
!   This routine smoothes an arbitrary field (fin) of length idim by applying
!   a digital filters of length nlength 4 nfilt times. The filterd field
!   is written on fout.
!
! Method:
!   Digital filter according to Shapiro
!
!------------------------------------------------------------------------------
!+ Implementation for double precision
!------------------------------------------------------------------------------

SUBROUTINE dfilt4_double (fin, idim, fhelp, fout, nfilt)

!------------------------------------------------------------------------------

! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)          ::    &
  idim,           & ! Dimension of the field
  nfilt             ! Number of iterative filerings
REAL (KIND=idouble),   INTENT (IN)          ::    &
  fin (idim)        ! input field (unfilterd)
REAL (KIND=idouble),   INTENT (OUT)         ::    &
  fout (idim)       ! smoothed output field (filtered)
REAL (KIND=idouble),   INTENT (INOUT)       ::    &
  fhelp(idim)       ! additional storage supplied by the calling routine

! Local variables
INTEGER (KIND=iintegers) ::    &
  i,m,            & ! loop indicees
  nf_o2             ! nfilt/2

REAL (KIND=idouble)      ::    & 
  fw(5)             ! filter weights

!------------------------------------------------------------------------------
  DATA fw / -0.00390625_idouble, 0.03125_idouble, -0.109375_idouble,     &
             0.21875_idouble,    0.7265625_idouble /


! begin subroutine dfilt4_double

  nf_o2 = (nfilt+1)/2

  fout (:) = fin(:)
  fhelp(:) = fin(:)

  DO i = 2, idim-1
    fhelp(i) =   0.15_idouble*fout (i-1) + 0.7_idouble*fout (i)    &
               + 0.15_idouble*fout (i+1)
  ENDDO
  DO i = 2, idim-1
    fout (i) =   0.15_idouble*fhelp(i-1) + 0.7_idouble*fhelp(i)    &
               + 0.15_idouble*fhelp(i+1)
  ENDDO

  DO m = 1, nf_o2
    DO i = 5, idim-4
      fhelp(i) =  fw(5)*fout(i) &
                + fw(4)*(fout(i-1)+fout(i+1)) + fw(3)*(fout(i-2)+fout(i+2)) &
                + fw(2)*(fout(i-3)+fout(i+3)) + fw(1)*(fout(i-4)+fout(i+4))  
    ENDDO
    DO i = 5, idim-4
      fout(i) = fw(5)*fhelp(i) &
              + fw(4)*(fhelp(i-1)+fhelp(i+1)) + fw(3)*(fhelp(i-2)+fhelp(i+2)) &
              + fw(2)*(fhelp(i-3)+fhelp(i+3)) + fw(1)*(fhelp(i-4)+fhelp(i+4))  
    ENDDO
  ENDDO

  DO i = 2, idim-1
    fhelp(i) =   0.15_idouble*fout (i-1) + 0.7_idouble*fout (i)   &
               + 0.15_idouble*fout (i+1)
  ENDDO
  DO i = 2, idim-1
    fout (i) =   0.15_idouble*fhelp(i-1) + 0.7_idouble*fhelp(i)   &
               + 0.15_idouble*fhelp(i+1)
  ENDDO

END SUBROUTINE dfilt4_double

!------------------------------------------------------------------------------
!+ Implementation for single precision
!------------------------------------------------------------------------------

SUBROUTINE dfilt4_single (fin, idim, fhelp, fout, nfilt)

!------------------------------------------------------------------------------

! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)          ::    &
  idim,           & ! Dimension of the field
  nfilt             ! Number of iterative filerings
REAL (KIND=isingle),   INTENT (IN)          ::    &
  fin (idim)        ! input field (unfilterd)
REAL (KIND=isingle),   INTENT (OUT)         ::    &
  fout (idim)       ! smoothed output field (filtered)
REAL (KIND=isingle),   INTENT (INOUT)       ::    &
  fhelp(idim)       ! additional storage supplied by the calling routine

! Local variables
INTEGER (KIND=iintegers) ::    &
  i,m,            & ! loop indicees
  nf_o2             ! nfilt/2

REAL (KIND=isingle)      ::    & 
  fw(5)             ! filter weights

!------------------------------------------------------------------------------
  DATA fw / -0.00390625_isingle, 0.03125_isingle, -0.109375_isingle,     &
             0.21875_isingle,    0.7265625_isingle /


! begin subroutine dfilt4_single

  nf_o2 = (nfilt+1)/2

  fout (:) = fin(:)
  fhelp(:) = fin(:)

  DO i = 2, idim-1
    fhelp(i) =   0.15_isingle*fout (i-1) + 0.7_isingle*fout (i)    &
               + 0.15_isingle*fout (i+1)
  ENDDO
  DO i = 2, idim-1
    fout (i) =   0.15_isingle*fhelp(i-1) + 0.7_isingle*fhelp(i)    &
               + 0.15_isingle*fhelp(i+1)
  ENDDO

  DO m = 1, nf_o2
    DO i = 5, idim-4
      fhelp(i) =  fw(5)*fout(i) &
                + fw(4)*(fout(i-1)+fout(i+1)) + fw(3)*(fout(i-2)+fout(i+2)) &
                + fw(2)*(fout(i-3)+fout(i+3)) + fw(1)*(fout(i-4)+fout(i+4))  
    ENDDO
    DO i = 5, idim-4
      fout(i) = fw(5)*fhelp(i) &
              + fw(4)*(fhelp(i-1)+fhelp(i+1)) + fw(3)*(fhelp(i-2)+fhelp(i+2)) &
              + fw(2)*(fhelp(i-3)+fhelp(i+3)) + fw(1)*(fhelp(i-4)+fhelp(i+4))  
    ENDDO
  ENDDO

  DO i = 2, idim-1
    fhelp(i) =   0.15_isingle*fout (i-1) + 0.7_isingle*fout (i)   &
               + 0.15_isingle*fout (i+1)
  ENDDO
  DO i = 2, idim-1
    fout (i) =   0.15_isingle*fhelp(i-1) + 0.7_isingle*fhelp(i)   &
               + 0.15_isingle*fhelp(i+1)
  ENDDO

END SUBROUTINE dfilt4_single

!------------------------------------------------------------------------------

!==============================================================================
!==============================================================================
!+ Defines all subroutines for the generic routine dfilt8
!------------------------------------------------------------------------------
!
! SUBROUTINE dfilt8 (fin, idim, fhelp, fout, nfilt)
!
!------------------------------------------------------------------------------
!
! Description:
!   This routine smoothes an arbitrary field (fin) of length idim by applying
!   a digital filters of length nlength 8 nfilt times. The filterd field
!   is written on fout.
!
! Method:
!   Digital filter according to Shapiro
!
!------------------------------------------------------------------------------
!+ Implementation for double precision
!------------------------------------------------------------------------------

SUBROUTINE dfilt8_double (fin, idim, fhelp, fout, nfilt)

!------------------------------------------------------------------------------

! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)          ::    &
  idim,           & ! Dimension of the field
  nfilt             ! Number of iterative filerings
REAL (KIND=idouble),   INTENT (IN)          ::    &
  fin (idim)        ! input field (unfilterd)
REAL (KIND=idouble),   INTENT (OUT)         ::    &
  fout (idim)       ! smoothed output field (filtered)
REAL (KIND=idouble),   INTENT (INOUT)       ::    &
  fhelp(idim)       ! additional storage supplied by the calling routine

! Local variables
INTEGER (KIND=iintegers) ::    &
  i,m,            & ! loop indicees
  nf_o2             ! nfilt/2

REAL (KIND=idouble)         ::  & 
  fw(9)  ! filter weights

!------------------------------------------------------------------------------
DATA fw /-0.0000152590_idouble,  0.0002441406_idouble, -0.0018310546_idouble, &
          0.0085449218_idouble, -0.0277709960_idouble,  0.0666503906_idouble, &
         -0.1221923828_idouble,  0.1745605469_idouble,  0.8036193848_idouble /

! begin subroutine dfilt8_double

  nf_o2 = (nfilt+1)/2

  fout (:) = fin(:)
  fhelp(:) = fin(:)

  DO i = 2, idim-1
    fhelp(i) =   0.25_idouble*fout (i-1) + 0.5_idouble*fout (i)     &
               + 0.25_idouble*fout (i+1)
  ENDDO
  DO i = 2, idim-1
    fout (i) =   0.25_idouble*fhelp(i-1) + 0.5_idouble*fhelp(i)     &
               + 0.25_idouble*fhelp(i+1)
  ENDDO

  DO m = 1, nf_o2
    DO i = 9, idim-8
      fhelp(i) = fw(9)*fout(i) &
               + fw(8)*(fout(i-1)+fout(i+1)) + fw(7)*(fout(i-2)+fout(i+2)) &
               + fw(6)*(fout(i-3)+fout(i+3)) + fw(5)*(fout(i-4)+fout(i+4)) &
               + fw(4)*(fout(i-5)+fout(i+5)) + fw(3)*(fout(i-6)+fout(i+6)) &
               + fw(2)*(fout(i-7)+fout(i+7)) + fw(1)*(fout(i-8)+fout(i+8))  
    ENDDO
    DO i = 9, idim-8
      fout(i) = fw(9)*fhelp(i) &
              + fw(8)*(fhelp(i-1)+fhelp(i+1)) + fw(7)*(fhelp(i-2)+fhelp(i+2)) &
              + fw(6)*(fhelp(i-3)+fhelp(i+3)) + fw(5)*(fhelp(i-4)+fhelp(i+4)) &
              + fw(4)*(fhelp(i-5)+fhelp(i+5)) + fw(3)*(fhelp(i-6)+fhelp(i+6)) &
              + fw(2)*(fhelp(i-7)+fhelp(i+7)) + fw(1)*(fhelp(i-8)+fhelp(i+8))
    ENDDO
  ENDDO

  DO i = 2, idim-1
    fhelp(i) =   0.25_idouble*fout (i-1) + 0.5_idouble*fout (i)    &
               + 0.25_idouble*fout (i+1)
  ENDDO
  DO i = 2, idim-1
    fout (i) =   0.25_idouble*fhelp(i-1) + 0.5_idouble*fhelp(i)    &
               + 0.25_idouble*fhelp(i+1)
  ENDDO

END SUBROUTINE dfilt8_double

!------------------------------------------------------------------------------
!+ Implementation for single precision
!------------------------------------------------------------------------------

SUBROUTINE dfilt8_single (fin, idim, fhelp, fout, nfilt)

!------------------------------------------------------------------------------

! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)          ::    &
  idim,           & ! Dimension of the field
  nfilt             ! Number of iterative filerings
REAL (KIND=isingle),   INTENT (IN)          ::    &
  fin (idim)        ! input field (unfilterd)
REAL (KIND=isingle),   INTENT (OUT)         ::    &
  fout (idim)       ! smoothed output field (filtered)
REAL (KIND=isingle),   INTENT (INOUT)       ::    &
  fhelp(idim)       ! additional storage supplied by the calling routine

! Local variables
INTEGER (KIND=iintegers) ::    &
  i,m,            & ! loop indicees
  nf_o2             ! nfilt/2

REAL (KIND=isingle)         ::  & 
  fw(9)  ! filter weights

!------------------------------------------------------------------------------
DATA fw /-0.0000152590_isingle,  0.0002441406_isingle, -0.0018310546_isingle, &
          0.0085449218_isingle, -0.0277709960_isingle,  0.0666503906_isingle, &
         -0.1221923828_isingle,  0.1745605469_isingle,  0.8036193848_isingle /

! begin subroutine dfilt8_single

  nf_o2 = (nfilt+1)/2

  fout (:) = fin(:)
  fhelp(:) = fin(:)

  DO i = 2, idim-1
    fhelp(i) =   0.25_isingle*fout (i-1) + 0.5_isingle*fout (i)     &
               + 0.25_isingle*fout (i+1)
  ENDDO
  DO i = 2, idim-1
    fout (i) =   0.25_isingle*fhelp(i-1) + 0.5_isingle*fhelp(i)     &
               + 0.25_isingle*fhelp(i+1)
  ENDDO

  DO m = 1, nf_o2
    DO i = 9, idim-8
      fhelp(i) = fw(9)*fout(i) &
               + fw(8)*(fout(i-1)+fout(i+1)) + fw(7)*(fout(i-2)+fout(i+2)) &
               + fw(6)*(fout(i-3)+fout(i+3)) + fw(5)*(fout(i-4)+fout(i+4)) &
               + fw(4)*(fout(i-5)+fout(i+5)) + fw(3)*(fout(i-6)+fout(i+6)) &
               + fw(2)*(fout(i-7)+fout(i+7)) + fw(1)*(fout(i-8)+fout(i+8))  
    ENDDO
    DO i = 9, idim-8
      fout(i) = fw(9)*fhelp(i) &
              + fw(8)*(fhelp(i-1)+fhelp(i+1)) + fw(7)*(fhelp(i-2)+fhelp(i+2)) &
              + fw(6)*(fhelp(i-3)+fhelp(i+3)) + fw(5)*(fhelp(i-4)+fhelp(i+4)) &
              + fw(4)*(fhelp(i-5)+fhelp(i+5)) + fw(3)*(fhelp(i-6)+fhelp(i+6)) &
              + fw(2)*(fhelp(i-7)+fhelp(i+7)) + fw(1)*(fhelp(i-8)+fhelp(i+8))
    ENDDO
  ENDDO

  DO i = 2, idim-1
    fhelp(i) =   0.25_isingle*fout (i-1) + 0.5_isingle*fout (i)    &
               + 0.25_isingle*fout (i+1)
  ENDDO
  DO i = 2, idim-1
    fout (i) =   0.25_isingle*fhelp(i-1) + 0.5_isingle*fhelp(i)    &
               + 0.25_isingle*fhelp(i+1)
  ENDDO

END SUBROUTINE dfilt8_single

!==============================================================================
!==============================================================================
!------------------------------------------------------------------------------

SUBROUTINE dolph (deltat, taus, m, window, t, time, time2, w, w2)

!------------------------------------------------------------------------------
!
! Description:
!  Calculation of Dolph-Chebyshev window or, for short, Dolph Window, using
!  the expression in the reference:
!    Antoniou, Andreas, 1993: Digital Filters: Analysis,
!    Design and Applications. McGraw-Hill, Inc., 689pp.
!
!  The Dolph window is optimal in the following sense:
!  For a given main-lobe width, the stop-band attenuation is minimal;
!  for a given stop-band level, the main-lobe width is minimal.
!
! Method:
!
! Modules used:    NONE
!
!------------------------------------------------------------------------------

! Parameter List:
! ---------------

INTEGER (KIND=iintegers), INTENT (IN)             ::  &
  m                   ! for dimensioning the work arrays

REAL  (KIND=ireals), INTENT (IN)                  ::  &
  deltat, taus        ! time step and cutoff period for filtering

REAL  (KIND=ireals), INTENT (OUT)                 ::  &
  window(0:2*m)       ! result

! The following variables are only used for work space
REAL  (KIND=ireals), INTENT (OUT)                 ::  &
  t(0:2*m), time(0:2*m), time2(0:2*m), w(0:2*m), w2(0:2*m)

! Local Variables:
! ----------------

INTEGER (KIND=iintegers)  :: nt, i, n, nm1, nn
REAL    (KIND=ireals)     :: zpi, zthetas, zx0, zarg, zterm1, zterm2, zrr,   &
                             zr, zdb, zsum, zsumw

!------------ End of header ---------------------------------------------------

! Begin subroutine dolph

  zpi = 4.0_ireals * ATAN(1.0_ireals)

  n = 2*m+1
  nm1 = n-1
  zthetas = 2.0_ireals*zpi*deltat/taus
  zx0 = 1.0_ireals / COS(zthetas/2.0_ireals)
  zterm1 = (zx0 + SQRT(zx0**2-1))**(REAL (N-1, ireals))
  zterm2 = (zx0 - SQRT(zx0**2-1))**(REAL (N-1, ireals))
  zrr = 0.5*(zterm1 + zterm2)
  zr = 1/zrr
  zdb = 20.0_ireals * LOG10(zr)

!------------------------------------------------------------

  DO nt = 0, M
    zsum = 1
    DO i = 1, M
      zarg = zx0 * cos(i*zpi/N)
      ! Calculate the Chebyshev polynomials
      ! Reference: Numerical Recipes, Page 184, recurrence
      !   T_n(x) = 2xT_{n-1}(x) - T_{n-2}(x) ,  n>=2.
      T(0) = 1
      T(1) = zarg
      DO nn=2,nm1
        T(nn) = 2*zarg*T(nn-1) - T(nn-2)
      ENDDO
      zterm1 = T(nm1)
      zterm2 = cos(2*nt*zpi*i/n)
      zsum   = zsum + zr*2 * zterm1 * zterm2
    ENDDO
    w(nt) = zsum / n
    TIME(nt) = nt
  ENDDO

  ! Fill in the negative-time values by symmetry.
  DO nt = 0, m
    w2(m+nt) = w(nt)
    w2(m-nt) = w(nt)
    time2(m+nt) =  time(nT)
    time2(m-nt) = -time(nT)
  ENDDO

  ! Fill up the array for return
  zsumw = 0.0_ireals
  DO nt = 0, 2*m
    zsumw = zsumw + w2(nt)
  ENDDO

  DO nt=0,2*m
    WINDOW(nt) = w2(nt)
  ENDDO
!
!
!----------------------------------------------------------
!       PRINT *, (w2(nT),    nT=0,2*M)
!----------------------------------------------------------
!

END SUBROUTINE dolph

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE elapsed_time    (realtimedif, istat)

!------------------------------------------------------------------------------
!
! Description:
!   Returns the elapsed wall-clock time in seconds since the last call. On
!   the first call the variables are only initialized. If no system clock is
!   present, an error value of istat=1 will be returned, if the optional
!   argument istat was passed from the calling routine. 
!   realtimedif is set to 0 then.
!
! Method:
!   The intrinsic function SYSTEM_CLOCK is used, that returns the number of
!   clock counts since some system dependent event in the past (e.g. midnight
!   for a 24-hour system clock). The difference of clock counts since the last
!   call is determined and converted into seconds. The variables "lfirst"
!   and "icountsold" (see below) have to be SAVEd for the next call.
!
! Modules used:    NONE
!
!------------------------------------------------------------------------------
!
! Parameter List:
! ---------------

REAL  (KIND=ireals), INTENT (OUT)                 ::  &
      realtimedif     ! wall-clock time since the last call in seconds
                      ! (0 if no system-clock is available)

INTEGER (KIND=iintegers), INTENT (OUT), OPTIONAL  ::  &
      istat           ! optional argument for error value


! Local Variables:
! ----------------

LOGICAL, SAVE      :: lfirst = .TRUE.   ! determine whether first call or not

INTEGER, SAVE      :: icountsold        ! number of counts in the last call

INTEGER            :: icountsnew,     & ! number of counts in this call
                      ir, im            ! other arguments to SYSTEM_CLOCK

LOGICAL            :: lpres             ! if optional argument is present

!------------ End of header ---------------------------------------------------

! Begin subroutine elapsed_time

  lpres = PRESENT (istat)

  CALL SYSTEM_CLOCK ( COUNT=icountsnew, COUNT_RATE=ir, COUNT_MAX=im )

  IF ( ir /= 0 ) THEN
    ! system clock is present
    IF (lpres) THEN
      istat = 0
    ENDIF

    IF (lfirst) THEN
      ! first call: store value for the number of clock counts
      icountsold = icountsnew
      lfirst     = .FALSE.
    ELSE
      ! convert the clock counts to seconds
      IF ( icountsnew >= icountsold ) THEN
        realtimedif = ( REAL (icountsnew - icountsold, ireals) )      &
                      / REAL (ir,ireals)
      ELSE
        realtimedif = REAL (im- (icountsold-icountsnew ), ireals)     &
                      / REAL (ir, ireals)
      ENDIF
      icountsold = icountsnew
    ENDIF
  ELSE
    ! no system clock present: set error value
    realtimedif = 0.0
    IF ( lpres ) THEN
      istat = 1
    ENDIF
  ENDIF

END SUBROUTINE elapsed_time

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE get_utc_date (ntsteps, ystartdate, dt, itype_calendar,           &
                         yactdate1, yactdate2, nactday, acthour)

!------------------------------------------------------------------------------
!
! Description:
!   This routine determines the actual date of this forecast step.
!
! Method:
!   Using the date of the forecast-start, the number of time steps 
!   already performed and the length of the time steps, the actual
!   date is calculated taking leap-years into consideration.
!   The date is given in three different formats.
!
! Modules used:    NONE
!
!------------------------------------------------------------------------------
!
! Input Parameter list:
! ---------------------

INTEGER   (KIND=iintegers), INTENT(IN)   ::                           &
  itype_calendar,   & ! for specifying the calendar used
  ntsteps             ! number of actual performed time-steps

REAL      (KIND=ireals), INTENT(IN)      ::                           &
  dt         ! time step in seconds

CHARACTER (LEN=10), INTENT(IN)           ::                           &
  ystartdate ! start date of the forecast

! Output Parameter list:
! ----------------------

CHARACTER (LEN=10), INTENT(OUT)          ::                           &
  yactdate1  ! actual date in the form   yyyymmddhh

CHARACTER (LEN=22), INTENT(OUT)          ::                           &
  yactdate2  ! actual date in the form   wd   dd.mm.yy  hh UTC


INTEGER   (KIND=iintegers), INTENT(OUT)  ::                           &
  nactday    ! day of the year

REAL      (KIND=ireals), INTENT(OUT)     ::                           &
  acthour    ! actual hour of the day

! Local variables:
INTEGER   (KIND=iintegers)   ::                                       &
  month(12), monthsum(13), ileap, iweek, iy, m,                       &
  idd, imm, iyy, ihh, iday, imonth, iyear, ihour, immhours, iyyhours, &
  iyear_hours

CHARACTER (LEN=3)            :: yweek(7)

! And for computing the amount of seconds of the whole forecast time,
! an 8-Byte INTEGER has to be used. Otherwise the computation fails after
! approx. 68 years!!

INTEGER, PARAMETER :: int_dp = KIND(1_8)
REAL(KIND=ireals)  :: zseconds

!------------ End of header ---------------------------------------------------

! Begin subroutine get_utc_date

DATA         month  / 31 ,  28 ,  31 ,  30 ,  31 ,  30 ,       &
                      31 ,  31 ,  30 ,  31 ,  30 ,  31 /
DATA         yweek  /'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN' /


! Statementfunction: ileap(yy) = 0:  no leap year, 
!                    ileap(yy) = 1:  leap year
! corrected version for Gregorian / Proleptic calendar
! by A. Dobler, CLM Community
  ileap (iy) = IABS( MOD(iy,  4) -   4) /   4  & ! every       4 years is a leapyear
              -IABS( MOD(iy,100) - 100) / 100  & ! every     100 years is no leapyear
              +IABS( MOD(iy,400) - 400) / 400    ! but every 400 years is a leapyear

! Divide ystartdate in day, month, year and hour
! and calculate the sums of days from the beginning of the year to the 
! end of the months
  READ ( ystartdate, '(I4,3I2)' ) iyy, imm, idd, ihh

  IF     (itype_calendar == 0) THEN
    month (2)    = 28 + ileap (iyy)
    monthsum(1) =  0
    DO m =  2 , 13
      monthsum(m) =  monthsum(m-1) + month(m-1)
    enddo
  ELSEIF (itype_calendar == 1) THEN
    monthsum(1) =  0
    DO m =  2 , 13
      monthsum(m) =  monthsum(m-1) + 30
    enddo
  ENDIF

! Determine how many hours have passed in this year
  iyyhours = (idd*24) + monthsum(imm)*24 + (ihh-24)
  iyyhours = iyyhours +                                           &
             INT (NINT (ntsteps*dt, int_dp)/3600_int_dp, iintegers)

! Take turning of the year into account
  IF     (itype_calendar == 0) THEN
    iyear_hours = 8760 + ileap(iyy)*24.0_ireals
  ELSEIF (itype_calendar == 1) THEN
    iyear_hours = 8640
  ENDIF

  IF (iyyhours < 0) THEN
    iyear    = iyy-1
    IF     (itype_calendar == 0) THEN
      iyyhours = 8760 + ileap(iyear)*24.0_ireals + iyyhours
    ELSEIF (itype_calendar == 1) THEN
      iyyhours = 8640                            + iyyhours
    ENDIF
  ELSE IF (iyyhours >= iyear_hours) THEN
    ! Take also into account if the run lasts for several years
    iyear    = iyy
    IF     (itype_calendar == 0) THEN
      iyear_hours = 8760 + ileap(iyear)*24.0_ireals
    ELSEIF (itype_calendar == 1) THEN
      iyear_hours = 8640
    ENDIF

    DO WHILE (iyyhours >= iyear_hours)
      iyyhours = iyyhours - iyear_hours
      iyear=iyear+1
      IF     (itype_calendar == 0) THEN
        iyear_hours = 8760 + ileap(iyear)*24.0_ireals
      ELSEIF (itype_calendar == 1) THEN
        iyear_hours = 8640
      ENDIF
    ENDDO
  ELSE
    iyear    =   iyy
  ENDIF

  ! calculate monthsum for actual year
  IF     (itype_calendar == 0) THEN
    month (2)    = 28 + ileap (iyear)
    monthsum(1) =  0
    DO m =  2 , 13
      monthsum(m) =  monthsum(m-1) + month(m-1)
    enddo
  ELSEIF (itype_calendar == 1) THEN
    monthsum(1) =  0
    DO m =  2 , 13
      monthsum(m) =  monthsum(m-1) + 30
    enddo
  ENDIF

! Determine the actual date from iyyhours
  m        = 1
  immhours = iyyhours
  DO WHILE (immhours >= 0)
    m        = m+1
    immhours = iyyhours - monthsum(m) * 24
  ENDDO
  imonth   = m-1

  immhours = iyyhours - monthsum(imonth)*24
  iday     = immhours/24 + 1
  ihour    = MOD(immhours,24)


!US: This was before, when 3600.0/dt was an integer value
!  acthour  = REAL (ihour, ireals) + dt/3600._ireals*  &
!                   MOD(ntsteps,INT(3600./dt+0.01)) + 0.0001_ireals

! this is the more accurate computation
  zseconds = REAL (ntsteps, ireals) * dt / 3600.0_ireals
  acthour  = REAL (ihour, ireals) +                          &
                  (zseconds - REAL(INT(zseconds, int_dp), ireals))

  ihour    = INT(acthour)
  nactday  = monthsum(imonth) + iday + INT(acthour/24. + 0.0001)
  iweek    = MOD(monthsum(imonth) + iday + (iyear-1901) + (iyear-1901)/4, 7)+1

  WRITE ( yactdate1(1:4) , '(I4.4)' ) iyear
  WRITE ( yactdate1(5:6) , '(I2.2)' ) imonth
  WRITE ( yactdate1(7:8) , '(I2.2)' ) iday
  WRITE ( yactdate1(9:10), '(I2.2)' ) ihour

  IF     (itype_calendar == 0) THEN
    yactdate2 = yweek(iweek)//' '//yactdate1(7:8)//'.'// yactdate1(5:6)//'.' &
                      //yactdate1(1:4)//'  '//yactdate1(9:10)//' UTC'
  ELSEIF (itype_calendar == 1) THEN
    yactdate2 = '    '//yactdate1(7:8)//'.'// yactdate1(5:6)//'.' &
                      //yactdate1(1:4)//'  '//yactdate1(9:10)//' UTC'
  ENDIF

END SUBROUTINE get_utc_date

!==============================================================================
!==============================================================================
!+ filter routines with 17-, 13-, 9- and 3-points stencil, resp.
!------------------------------------------------------------------------------

SUBROUTINE horizontal_filtering( field_flt, ie_in, je_in, kedim,          &
                                 nbdlines, nflt_width, ncutoff,           &
                                 neighbors, hfx_mask, hfy_mask )

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine arguments:
! ---------------------
INTEGER (KIND=iintegers), INTENT(IN) :: &
  ie_in, je_in, kedim,  & ! Dimensions of the field to be filtered
  nbdlines,             & ! number of boundary lines from decomposition
  nflt_width,           & ! width of field extension for filtering
  ncutoff,              & ! filter-value for cutoff
  neighbors(4)            ! process-id's of the neighbors in the grid


REAL    (KIND=ireals   ), INTENT(INOUT) ::  &
  field_flt(ie_in, je_in, kedim)

LOGICAL, INTENT(in), OPTIONAL ::  &
  hfx_mask(ie_in, je_in), hfy_mask(ie_in, je_in)

! Local scalars:
! -------------
INTEGER (KIND=iintegers) ::  &
  ilow, iup,           & !
  jlow, jup,           & !
  izstata,             & !  error status at allocation
  izstatd,             & !  error status at deallocation
  i, j, k, l             !  Loop indices

INTEGER (KIND=iintegers) ::  &
  istart, iend, jstart, jend, nfw_m_nb

! Local (automatic) arrays:
! -------------------------
REAL    (KIND=ireals   ) ::  &
  field_tmp (ie_in, je_in, kedim), &
  field_tmp2(ie_in, je_in, kedim), &
  zfwnp(-nflt_width:nflt_width),   & ! filter weights for n-point filter
  zfw3p(-1:1)                        ! filter weights for 3-point filter

!------------------------------------------------------------------------------

  nfw_m_nb = nflt_width - nbdlines
  istart = 1 + nbdlines
  iend   = ie_in - 2*nfw_m_nb - nbdlines
  jstart = 1 + nbdlines
  jend   = je_in - 2*nfw_m_nb - nbdlines

  ! filter weights for n-point filter
  IF (ncutoff == 3 .AND. nflt_width == 4) THEN
    ! --> dfilt4
    ! filter weights for 9-point filter (approx. cutoff = 3)
    zfwnp = (/ -0.390625E-02_ireals,     &
               +0.3125E-01_ireals,       &
               -0.109375_ireals,         &
               +0.21875_ireals,          &
               +0.7265625_ireals,        &
               +0.21875_ireals,          &
               -0.109375_ireals,         &
               +0.3125E-01_ireals,       &
               -0.390625E-02_ireals /)
  ELSEIF (ncutoff == 3 .AND. nflt_width == 8) THEN
    ! --> dfilt8
    ! filter weights for 17-point filter (approx. cutoff = 3)
    zfwnp = (/ -0.15259E-04_ireals,      &
               +0.2441406E-03_ireals,    &
               -0.18310546E-02_ireals,   &
               +0.85449218E-02_ireals,   &
               -0.27770996E-01_ireals,   &
               +0.666503906E-01_ireals,  &
               -0.1221923828_ireals,     &
               +0.1745605469_ireals,     &
               +0.8036193848_ireals,     &
               +0.1745605469_ireals,     &
               -0.1221923828_ireals,     &
               +0.666503906E-01_ireals,  &
               -0.27770996E-01_ireals,   &
               +0.85449218E-02_ireals,   &
               -0.18310546E-02_ireals,   &
               +0.2441406E-03_ireals,    &
               -0.15259E-04_ireals /)
  ELSEIF (ncutoff == 4 .AND. nflt_width == 4) THEN
    ! filter weights for 9-point filter (approx. cutoff = 4)
    zfwnp = (/ +0.1171875E-01_ireals,    &
               -0.3125E-01_ireals,       &
               -0.46875E-01_ireals,      &
               +0.28125_ireals,          &
               +0.5703125_ireals,        &
               +0.28125_ireals,          &
               -0.46875E-01_ireals,      &
               -0.3125E-01_ireals,       &
               +0.1171875E-01_ireals /)
  ELSEIF (ncutoff == 5 .AND. nflt_width == 6) THEN
    ! filter weights for 13-point filter (approx. cutoff = 5)
    zfwnp = (/ +0.44023278E-02_ireals,   &
               +0.13175894E-01_ireals,   &
               -0.477203075E-01_ireals,  &
               -0.435555245E-01_ireals,  &
               +0.94700467E-01_ireals,   &
               +0.2888298641_ireals,     &
               +0.3803345582_ireals,     &
               +0.2888298641_ireals,     &
               +0.94700467E-01_ireals,   &
               -0.435555245E-01_ireals,  &
               -0.477203075E-01_ireals,  &
               +0.13175894E-01_ireals,   &
               +0.44023278E-02_ireals /)
  ELSEIF (ncutoff == 6 .AND. nflt_width == 4) THEN
    ! filter weights for 9-point filter (approx. cutoff = 6)
    zfwnp = (/ -0.4694126E-01_ireals,    &
               -0.50095541E-02_ireals,   &
               +0.13528415_ireals,       &
               +0.25500955_ireals,       &
               +0.32331423_ireals,       &
               +0.25500955_ireals,       &
               +0.13528415_ireals,       &
               -0.50095541E-02_ireals,   &
               -0.4694126E-01_ireals /)
  ELSEIF (ncutoff == 8 .AND. nflt_width == 6) THEN
    ! filter weights for 13-point filter (approx. cutoff = 8)
    zfwnp = (/ -0.16638111E-01_ireals,   &
               -0.30753028E-01_ireals,   &
               -0.17361869E-02_ireals,   &
               +0.65428931E-01_ireals,   &
               +0.14784805_ireals,       &
               +0.2153241_ireals,        &
               +0.2410525_ireals,        &
               +0.2153241_ireals,        &
               +0.14784805_ireals,       &
               +0.65428931E-01_ireals,   &
               -0.17361869E-02_ireals,   &
               -0.30753028E-01_ireals,   &
               -0.16638111E-01_ireals /)
  ELSE
    PRINT *, ' ERROR *** Wrong cutoff value for filtering        or *** '
    PRINT *, ' ERROR *** wrong value for filter/field extension.    *** '
  ENDIF

  ! filter weights for 3-point filter (approx. cutoff = 4)
  zfw3p = (/ 0.25_ireals, 0.5_ireals, 0.25_ireals /)

  ! west
  IF (neighbors(1) == -1) THEN
    ilow = 1 + 2*nflt_width
  ELSE
    ilow = istart + nfw_m_nb
  END IF
  ! east
  IF (neighbors(3) == -1) THEN
    iup = iend - nbdlines
  ELSE
    iup = iend + nfw_m_nb
  END IF
  ! south
  IF (neighbors(4) == -1) THEN
    jlow = 1 + 2*nflt_width
  ELSE
    jlow = jstart + nfw_m_nb
  END IF
  ! north
  IF (neighbors(2) == -1) THEN
    jup = jend - nbdlines
  ELSE
    jup = jend + nfw_m_nb
  END IF

  ! init working array
  field_tmp (:,:,:) = field_flt(:,:,:)


  IF ( PRESENT( hfx_mask ) ) THEN

    ! apply n-point-filter in x-direction
    DO k = 1, kedim
      DO j = 1, je_in
        DO i = ilow, iup
          IF ( hfx_mask(i,j) ) THEN
            field_tmp(i,j,k) = 0.0_ireals
          ENDIF
        ENDDO
        DO l = -nflt_width, nflt_width
          DO i = ilow, iup
            IF ( hfx_mask(i,j) ) THEN
              field_tmp(i,j,k) = field_tmp(i,j,k)               &
                               + zfwnp(l)*field_flt(i+l,j,k)
            END IF
          END DO
        END DO
      END DO
    END DO

    ! apply 3-point-filter in x-direction at west boundary
    IF (neighbors(1) == -1) THEN
      DO k = 1, kedim
        DO j = 1, je_in
          DO i = nfw_m_nb+1, ilow-1
            IF ( hfx_mask(i,j) ) THEN
              field_tmp(i,j,k) = 0.0_ireals
            ENDIF
          ENDDO
          DO l = -1, 1
            DO i = nfw_m_nb+1, ilow-1
              IF ( hfx_mask(i,j) ) THEN
                field_tmp(i,j,k) = field_tmp(i,j,k)             &
                                 + zfw3p(l)*field_flt(i+l,j,k)
              END IF
            END DO
          END DO
        END DO
      END DO
    END IF

    ! apply 3-point-filter in x-direction at east boundary
    IF (neighbors(3) == -1) THEN
      DO k = 1, kedim
        DO j = 1, je_in
          DO i = iup+1, ie_in-nfw_m_nb
            IF ( hfx_mask(i,j) ) THEN
              field_tmp(i,j,k) = 0.0_ireals
            ENDIF
          ENDDO
          DO l = -1, 1
            DO i = iup+1, ie_in-nfw_m_nb
              IF ( hfx_mask(i,j) ) THEN
                field_tmp(i,j,k) = field_tmp(i,j,k)             &
                                 + zfw3p(l)*field_flt(i+l,j,k)
              END IF
            END DO
          END DO
        END DO
      END DO
    END IF

  ELSE

    !
    ! apply n-point-filter in x-direction
    !
    DO k = 1, kedim
      DO j = 1, je_in
        DO i = ilow, iup
          field_tmp(i,j,k) = 0.0_ireals
        END DO
        DO l = -nflt_width, nflt_width
          DO i = ilow, iup
            field_tmp(i,j,k) = field_tmp(i,j,k)                 &
                             + zfwnp(l)*field_flt(i+l,j,k)
          END DO
        END DO
      END DO
    END DO

    ! apply 3-point-filter in x-direction at west boundary
    IF (neighbors(1) == -1) THEN
      DO k = 1, kedim
        DO j = 1, je_in
          DO i = nfw_m_nb+1, ilow-1
            field_tmp(i,j,k) = 0.0_ireals
          END DO
          DO l = -1, 1
            DO i = nfw_m_nb+1, ilow-1
              field_tmp(i,j,k) = field_tmp(i,j,k)               &
                               + zfw3p(l)*field_flt(i+l,j,k)
            END DO
          END DO
        END DO
      END DO
    END IF

    ! apply 3-point-filter in x-direction at east boundary
    IF (neighbors(3) == -1) THEN
      DO k = 1, kedim
        DO j = 1, je_in
          DO i = iup+1, ie_in-nfw_m_nb
            field_tmp(i,j,k) = 0.0_ireals
          END DO
          DO l = -1, 1
            DO i = iup+1, ie_in-nfw_m_nb
              field_tmp(i,j,k) = field_tmp(i,j,k)               &
                               + zfw3p(l)*field_flt(i+l,j,k)
            END DO
          END DO
        END DO
      END DO
    END IF

  END IF


  IF ( PRESENT( hfy_mask ) ) THEN

    ! apply n-point-filter in y-direction
    DO k = 1, kedim
      DO j = jlow, jup
        DO i = 1, ie_in
          IF ( hfy_mask(i,j) ) THEN
            field_flt(i,j,k) = 0.0_ireals
          ELSE
            field_flt(i,j,k) = field_tmp(i,j,k)
          ENDIF
        ENDDO
        DO l = -nflt_width, nflt_width
          DO i = 1, ie_in
            IF ( hfy_mask(i,j) ) THEN
              field_flt(i,j,k) = field_flt(i,j,k) + zfwnp(l)*field_tmp(i,j+l,k)
            END IF
          END DO
        END DO
      END DO
    END DO

    ! apply 3-point-filter in y-direction at south boundary
    IF (neighbors(4) == -1) THEN
      DO k = 1, kedim
        DO j = nfw_m_nb+1, jlow-1
          DO i = 1, ie_in
            IF ( hfy_mask(i,j) ) THEN
              field_flt(i,j,k) = 0.0_ireals
            ELSE
              field_flt(i,j,k) = field_tmp(i,j,k)
            ENDIF
          ENDDO
          DO l = -1, 1
            DO i = 1, ie_in
              IF ( hfy_mask(i,j) ) THEN
                field_flt(i,j,k) = field_flt(i,j,k)+zfw3p(l)*field_tmp(i,j+l,k)
              END IF
            END DO
          END DO
        END DO
      END DO
    END IF

    ! apply 3-point-filter in y-direction at north boundary
    IF (neighbors(2) == -1) THEN
      DO k = 1, kedim
        DO j = jup+1, je_in-nfw_m_nb
          DO i = 1, ie_in
            IF ( hfy_mask(i,j) ) THEN
              field_flt(i,j,k) = 0.0_ireals
            ELSE
              field_flt(i,j,k) = field_tmp(i,j,k)
            ENDIF
          ENDDO
          DO l = -1, 1
            DO i = 1, ie_in
              IF ( hfy_mask(i,j) ) THEN
                field_flt(i,j,k) = field_flt(i,j,k)+zfw3p(l)*field_tmp(i,j+l,k)
              END IF
            END DO
          END DO
        END DO
      END DO
    END IF

  ELSE

    !
    ! apply n-point-filter in y-direction
    !
    DO k = 1, kedim
      DO j = jlow, jup
        DO i = 1, ie_in
          field_flt(i,j,k) = 0.0_ireals
        ENDDO
        DO l = -nflt_width, nflt_width
          DO i = 1, ie_in
            field_flt(i,j,k) = field_flt(i,j,k)+zfwnp(l)*field_tmp(i,j+l,k)
          END DO
        END DO
      END DO
    END DO

    ! apply 3-point-filter in y-direction at south boundary
    IF (neighbors(4) == -1) THEN
      DO k = 1, kedim
        DO j = nfw_m_nb+1, jlow-1
          DO i = 1, ie_in
            field_flt(i,j,k) = 0.0_ireals
          ENDDO
          DO l = -1, 1
            DO i = 1, ie_in
              field_flt(i,j,k) = field_flt(i,j,k)+zfw3p(l)*field_tmp(i,j+l,k)
            END DO
          END DO
        END DO
      END DO
    END IF

    ! apply 3-point-filter in y-direction at north boundary
    IF (neighbors(2) == -1) THEN
      DO k = 1, kedim
        DO j = jup+1, je_in-nfw_m_nb
          DO i = 1, ie_in
            field_flt(i,j,k) = 0.0_ireals
          ENDDO
          DO l = -1, 1
            DO i = 1, ie_in
              field_flt(i,j,k) = field_flt(i,j,k)+zfw3p(l)*field_tmp(i,j+l,k)
            END DO
          END DO
        END DO
      END DO
    END IF

  END IF

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! End of subroutine horizontal_filtering
!------------------------------------------------------------------------------

END SUBROUTINE horizontal_filtering

!==============================================================================
!==============================================================================
!+ Function for rotation of geographical coordinates
!------------------------------------------------------------------------------

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
REAL (KIND=ireals), INTENT (IN)      ::        &
  polphi,   & ! latitude of the rotated north pole
  pollam,   & ! longitude of the rotated north pole
  phirot,   & ! latitude in the rotated system
  rlarot      ! longitude in the rotated system

REAL (KIND=ireals), INTENT (IN)      ::        &
  polgam      ! angle between the north poles of the systems

REAL (KIND=ireals)                   ::        &
  phirot2phi  ! latitude in the geographical system

! Local variables
REAL (KIND=ireals)                   ::        &
  zsinpol, zcospol, zphis, zrlas, zarg, zgam

REAL (KIND=ireals), PARAMETER        ::        &
  zrpi18 = 57.2957795_ireals,                  &
  zpir18 = 0.0174532925_ireals

!------------------------------------------------------------------------------

! Begin function phirot2phi

  zsinpol     = SIN (zpir18 * polphi)
  zcospol     = COS (zpir18 * polphi)
 
  zphis       = zpir18 * phirot
  IF (rlarot > 180.0_ireals) THEN
    zrlas = rlarot - 360.0_ireals
  ELSE
    zrlas = rlarot
  ENDIF
  zrlas       = zpir18 * zrlas

  IF (polgam /= 0.0_ireals) THEN
    zgam  = zpir18 * polgam
    zarg  = zsinpol*SIN (zphis) +                                           &
        zcospol*COS(zphis) * ( COS(zrlas)*COS(zgam) - SIN(zgam)*SIN(zrlas) )
  ELSE
    zarg  = zcospol * COS (zphis) * COS (zrlas) + zsinpol * SIN (zphis)
  ENDIF
 
  phirot2phi  = zrpi18 * ASIN (zarg)

END FUNCTION phirot2phi

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

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
REAL (KIND=ireals), INTENT (IN)      ::        &
  polphi,  & ! latitude of the rotated north pole
  pollam,  & ! longitude of the rotated north pole
  phi,     & ! latitude in the geographical system
  rla        ! longitude in the geographical system

REAL (KIND=ireals)                   ::        &
  phi2phirot ! longitude in the rotated system

! Local variables
REAL (KIND=ireals)                       ::    &
  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

REAL (KIND=ireals), PARAMETER            ::    &
  zrpi18 = 57.2957795_ireals,                  & !
  zpir18 = 0.0174532925_ireals

!------------------------------------------------------------------------------

! Begin function phi2phirot

  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0_ireals) THEN
    zrla1  = rla - 360.0_ireals
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
REAL (KIND=ireals), INTENT (IN)      ::        &
  polphi,   & ! latitude of the rotated north pole
  pollam,   & ! longitude of the rotated north pole
  phirot,   & ! latitude in the rotated system
  rlarot      ! longitude in the rotated system

REAL (KIND=ireals), INTENT (IN)      ::        &
  polgam      ! angle between the north poles of the systems

REAL (KIND=ireals)                   ::        &
  rlarot2rla  ! latitude in the geographical system

! Local variables
REAL (KIND=ireals)                   ::        &
  zsinpol, zcospol, zlampol, zphis, zrlas, zarg1, zarg2, zgam

REAL (KIND=ireals), PARAMETER        ::        &
  zrpi18 = 57.2957795_ireals,                  & !
  zpir18 = 0.0174532925_ireals

!------------------------------------------------------------------------------

! Begin function rlarot2rla

  zsinpol = SIN (zpir18 * polphi)
  zcospol = COS (zpir18 * polphi)

  zlampol = zpir18 * pollam
  zphis   = zpir18 * phirot
  IF (rlarot > 180.0_ireals) THEN
    zrlas = rlarot - 360.0_ireals
  ELSE
    zrlas = rlarot
  ENDIF
  zrlas   = zpir18 * zrlas

  IF (polgam /= 0.0_ireals) THEN
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
 
  IF (zarg2 == 0.0) zarg2 = 1.0E-20_ireals
 
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
REAL (KIND=ireals), INTENT (IN)      ::        &
  polphi,  & ! latitude of the rotated north pole
  pollam,  & ! longitude of the rotated north pole
  phi,     & ! latitude in geographical system
  rla        ! longitude in geographical system

REAL (KIND=ireals), INTENT (IN)      ::        &
  polgam      ! angle between the north poles of the systems

REAL (KIND=ireals)                   ::        &
  rla2rlarot ! latitude in the the rotated system

! Local variables
REAL (KIND=ireals)                       ::    &
  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

REAL (KIND=ireals), PARAMETER            ::    &
  zrpi18 = 57.2957795_ireals,                  & !
  zpir18 = 0.0174532925_ireals

!------------------------------------------------------------------------------

! Begin function rla2rlarot

  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0_ireals) THEN
    zrla1  = rla - 360.0_ireals
  ELSE
    zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1

  zarg1    = - SIN (zrla-zlampol) * COS(zphi)
  zarg2    = - zsinpol * COS(zphi) * COS(zrla-zlampol) + zcospol * SIN(zphi)

  IF (zarg2 == 0.0) zarg2 = 1.0E-20_ireals

  rla2rlarot = zrpi18 * ATAN2 (zarg1,zarg2)

  IF (polgam /= 0.0_ireals ) THEN
    rla2rlarot = polgam + rla2rlarot
    IF (rla2rlarot > 180._ireals) rla2rlarot = rla2rlarot -360._ireals
  ENDIF

END FUNCTION rla2rlarot

!==============================================================================
!==============================================================================

SUBROUTINE sleve_split_oro (hsurf, hsurfs, idim, jdim, nflt, nextralines,   &
                            svc1, svc2, vcflat, noutunit, myid, ierror, yerror)

!------------------------------------------------------------------------------
!
! Description:
!     decomposes a given topography field hsurf in a
!     large-scale (hsurfs(:,:,1)) and a small-scale (hsurfs(:,:,2)) part, where
!     hsurf(:,:) = hsurfs(:,:,1) + hsurfs(:,:,2)
!
! Method:
!     - a digital filter is applied for the computation of
!       the large scale part hsurfs(:,:,1).
!     - the boundary values are treated seperately to assure, that
!       also these points are smoothed:
!       i.e. at the i=1    boundary: A(1,j)    = A(2,j)      for all j
!                   i=idim boundary: A(idim,j) = A(idim-1,j) for all j
!                   j=1    boundary: A(i,1)    = A(i,2)      for all i
!                   j=jdim boundary: A(i,jdim) = A(i,jdim-1) for all i
!     - nflt determines, how often the filter is applied
!     - Additionally, the maxima of hsurf, hsurfs(:,:,1) and hsurfs(:,:,2) are
!       computed and written to noutunit.
!
!    written by Daniel Leuenberger, 03.10.2001
!------------------------------------------------------------------------------

! Subroutine Arguments:

INTEGER  (KIND=iintegers), INTENT(IN)       ::    &
          idim, jdim,                & !  dimensions of hsurf
          nextralines,               & !  number of extra lines around filtered
                                       !  field (for interpolation program)
          nflt                         !  number of filter applications

REAL     (KIND=ireals), INTENT(IN)          ::    &
          svc1, svc2,                & !  decay rates for large and small scale
          vcflat                       !  vertical coordinate where the
                                       !  terrain following system changes back
                                       !  to an orthogonal z-system
REAL     (KIND=ireals), INTENT(IN)          ::    &
          hsurf(idim,jdim)             !  height of full topography

INTEGER  (KIND=iintegers), INTENT(IN)       ::    &
          noutunit                     !  unitnumber where output is written to

REAL     (KIND=ireals), INTENT(OUT)         ::    &
          hsurfs(idim,jdim,2)          !  height of splitted topography parts

INTEGER  (KIND=iintegers), INTENT(IN)       ::    &
          myid                         !  PE number

INTEGER  (KIND=iintegers), INTENT(OUT)      ::    &
          ierror                       !  error value

CHARACTER(LEN=*)                            ::    &
          yerror                       !  error message

! Local variables
REAL     (KIND=ireals)                      ::    &
          maxhsurf,                 &  !  maximum of hsurf
          maxhsurf1,                &  !  maximum of hsurfs(:,:,1)
          maxhsurf2,                &  !  maximum of hsurfs(:,:,2)
          gammavc                      !  invertibility parameter

INTEGER  (KIND=iintegers)                   ::    &
          i,j,n,old,new,temp, istart, iend, jstart, jend

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE sleve_split_oro
!------------------------------------------------------------------------------

  ierror    = 0_iintegers
  yerror    = '        '

  IF (myid == 0) THEN
    WRITE (noutunit,'(A)') '    '
    WRITE (noutunit,'(A)') '   Splitting of Topography for SLEVE coordinate:'
  ENDIF

  IF ( (nextralines < 0) .OR. (nextralines > 1) ) THEN
    ierror = 1
    yerror = 'ERROR:  nextralines outside range: 0 <= nextralines <= 1'
    RETURN
  ENDIF

  ! In order to obtain the same splitting in LM and in INT2LM, only the domain
  ! without the extra boundary lines is splitted. These extra boundary lines
  ! are indicated by nextralines

  ! compute index boundaries of LM-domain
  istart    = 1    + nextralines
  iend      = idim - nextralines
  jstart    = 1    + nextralines
  jend      = jdim - nextralines

  maxhsurf  = 0.0_ireals
  maxhsurf1 = 0.0_ireals
  maxhsurf2 = 0.0_ireals

  old = 1
  new = 2

  hsurfs(:,:,old) = hsurf(:,:)

  ! apply nflt times an ideal 2d filter to compute the
  ! large-scale part hsurfs(:,:,1) of the topography

  DO n = 1, nflt
    ! treat inner points
    DO j=jstart+1, jend-1
      DO i=istart+1, iend-1
        hsurfs(i,j,new) = 0.25_ireals * hsurfs(i,j,old)                      &
              + 0.125_ireals  * (hsurfs(i-1,j  ,old) + hsurfs(i+1,j  ,old) + &
                                 hsurfs(i  ,j-1,old) + hsurfs(i  ,j+1,old))  &
              + 0.0625_ireals * (hsurfs(i-1,j-1,old) + hsurfs(i+1,j-1,old) + &
                                 hsurfs(i-1,j+1,old) + hsurfs(i+1,j+1,old))
      ENDDO
    ENDDO

    ! treat corner points
    hsurfs(istart,jstart,new) =  0.25_ireals *                               &
            (hsurfs(istart, jstart  ,old) + hsurfs(istart+1, jstart  ,old)   &
           + hsurfs(istart, jstart+1,old) + hsurfs(istart+1, jstart+1,old))

    hsurfs(istart, jend, new) =  0.25_ireals *                               &
                 (hsurfs(istart, jend  ,old) + hsurfs(istart+1, jend  ,old)  &
                + hsurfs(istart, jend-1,old) + hsurfs(istart+1, jend-1,old))

    hsurfs(iend, jstart, new) =  0.25_ireals *                               &
                (hsurfs(iend, jstart  ,old) + hsurfs(iend-1,jstart  ,old)    &
               + hsurfs(iend, jstart+1,old) + hsurfs(iend-1,jstart+1,old))

    hsurfs(iend,jend,new) =  0.25_ireals *                                   &
                         (hsurfs(iend,jend  ,old) + hsurfs(iend-1,jend  ,old)&
                        + hsurfs(iend,jend-1,old) + hsurfs(iend-1,jend-1,old))

    ! treat edge points
    DO j = jstart+1,jend-1
      hsurfs(istart,j,new)  =                                                &
        0.25_ireals  * (hsurfs(istart  ,j  ,old) + hsurfs(istart+1,j  ,old)) &
     +  0.125_ireals * (hsurfs(istart  ,j-1,old) + hsurfs(istart  ,j+1,old) +&
                        hsurfs(istart+1,j-1,old) + hsurfs(istart+1,j+1,old) )


      hsurfs(iend,j,new) =                                                   &
         0.25_ireals  * (hsurfs(iend  ,j  ,old) + hsurfs(iend-1,j  ,old))    &
      +  0.125_ireals * (hsurfs(iend  ,j-1,old) + hsurfs(iend  ,j+1,old) +   &
                         hsurfs(iend-1,j-1,old) + hsurfs(iend-1,j+1,old) )
    ENDDO

    DO i = istart+1,iend-1
      hsurfs(i,jstart,new)  =                                                &
        0.25_ireals  * (hsurfs(i  ,jstart  ,old) + hsurfs(i  ,jstart+1,old)) &
     +  0.125_ireals * (hsurfs(i-1,jstart  ,old) + hsurfs(i+1,jstart  ,old) +&
                        hsurfs(i-1,jstart+1,old) + hsurfs(i+1,jstart+1,old) )

      hsurfs(i,jend,new) =                                                   &
         0.25_ireals  * (hsurfs(i  ,jend  ,old) + hsurfs(i  ,jend-1,old))    &
      +  0.125_ireals * (hsurfs(i-1,jend  ,old) + hsurfs(i+1,jend  ,old) +   &
                         hsurfs(i-1,jend-1,old) + hsurfs(i+1,jend-1,old) )
    ENDDO

    temp = old
    old  = new
    new  = temp

  ENDDO

  ! compute the large-scale part hsurfs(:,:,1) of the topo
  hsurfs(istart:iend,jstart:jend,1) = hsurfs(istart:iend,jstart:jend,old)

  ! compute the small-scale part hsurfs(:,:,2) of the topo
  hsurfs(istart:iend,jstart:jend,2) = hsurf (istart:iend,jstart:jend) -    &
                                      hsurfs(istart:iend,jstart:jend,1)

  ! compute maxima of topographies
  maxhsurf  = MAXVAL (hsurf (istart:iend,jstart:jend)  )
  maxhsurf1 = MAXVAL (hsurfs(istart:iend,jstart:jend,1))
  maxhsurf2 = MAXVAL (hsurfs(istart:iend,jstart:jend,2))

  IF (myid == 0) THEN
    WRITE(noutunit,'(A,I5,A)' ) '    nflt = ',nflt,' Applications of Filter'
    WRITE(noutunit,'(A)'      ) '    Maxima of Topography Parts:'
    WRITE(noutunit,'(A,F10.5)')                                               &
             '    Max of Full Topography        hsurf          : ',maxhsurf
    WRITE(noutunit,'(A,F10.5)')                                               &
             '    Max of Large-Scale Topography hsurfs(:,:,1)  : ',maxhsurf1
    WRITE(noutunit,'(A,F10.5)')                                               &
             '    Max of Small-Scale Topography hsurfs(:,:,2)  : ',maxhsurf2
  ENDIF

  ! calculate SLEVE invertibility parameter gammavc
  gammavc = 1.0_ireals -                                                      &
    ((MAXVAL(hsurfs(istart:iend,jstart:jend,1)) / svc1) / TANH(vcflat/svc1))- &
    ((MAXVAL(hsurfs(istart:iend,jstart:jend,2)) / svc2) / TANH(vcflat/svc2))

  IF (myid == 0) THEN
    WRITE (noutunit,'(A)') '         '
    WRITE (noutunit,'(A,F10.5)')                                      &
       '   Invertibility parameter for SLEVE coordinate: gammavc = ',gammavc
    WRITE (noutunit,'(A)') '         '
  ENDIF

  ! check if invertibility condition is fulfilled
  IF ( gammavc <= 0.0 ) THEN
    PRINT *, 'Invertibility parameter for SLEVE coordinate: gammavc = ',&
              gammavc
    PRINT *, 'vcflat = ',vcflat
    PRINT *, 'svc1   = ',svc1
    PRINT *, 'svc2   = ',svc2
    ierror  = 2
    yerror  = 'Invertibility condition of SLEVE coordinate not '// &
              'fulfilled, check values of svc1, svc2 and vcflat'
    RETURN
  ELSEIF ( gammavc < 0.05 ) THEN
    PRINT *, 'Invertibility parameter for SLEVE coordinate: gammavc = ',&
              gammavc
    PRINT *, 'vcflat = ',vcflat
    PRINT *, 'svc1   = ',svc1
    PRINT *, 'svc2   = ',svc2
    PRINT *, 'WARNING !!! SLEVE Invertibility parameter close to ',   &
             'zero, check values of svc1, svc2 and vcflat'
  ENDIF

  IF (nextralines > 0) THEN
    ! The values of hsurfs outside the LM-domain are determined as follows:
    ! First, the large-scale topo hsurfs(:,:,1) is linearly extrapolated
    !     (from the point at the boundary of the LM-domain and the
    !      first point inside the LM-domain),
    ! then the small-scale topo hsurfs(:,:,2) is calculated from the
    ! relationship h2 = h - h1

    ! ** Attention:  This extrapolation works only in the case of
    ! ** nextralines = 1 !!!
    ! ** For nextralines > 1 the extrapolation has yet to be implemented !!!

    DO i = istart, iend
      ! extrapolation of south edge
      hsurfs(i,1,1) = 2 * hsurfs(i,2,1) - hsurfs(i,3,1)
      hsurfs(i,1,2) =     hsurf (i,1)   - hsurfs(i,1,1)

      ! extrapolation of north edge
      hsurfs(i,jdim,1) = 2 * hsurfs(i,jdim-1,1) - hsurfs(i,jdim-2,1)
      hsurfs(i,jdim,2) =     hsurf (i,jdim)     - hsurfs(i,jdim  ,1)
    ENDDO

    DO j = jstart, jend
      ! extrapolation of west edge
      hsurfs(1,j,1) = 2 * hsurfs(2,j,1) - hsurfs(3,j,1)
      hsurfs(1,j,2) =     hsurf (1,j)   - hsurfs(1,j,1)

      ! extrapolation of east edge
      hsurfs(idim,j,1) = 2 * hsurfs(idim-1,j,1) - hsurfs(idim-2,j,1)
      hsurfs(idim,j,2) =     hsurf (idim  ,j)   - hsurfs(idim  ,j,1)
    ENDDO

    ! extrapolation of SW point
    hsurfs(1,1,1) = 2 * hsurfs(2,2,1) - hsurfs(3,3,1)
    hsurfs(1,1,2) = hsurf(1,1) - hsurfs(1,1,1)

    ! extrapolation of  SE point
    hsurfs(idim,1,1) = 2 * hsurfs(idim-1,2,1) - hsurfs(idim-2,3,1)
    hsurfs(idim,1,2) =     hsurf (idim  ,1)   - hsurfs(idim  ,1,1)

    ! extrapolation of NW point
    hsurfs(1,jdim,1) = 2 * hsurfs(2,jdim-1,1) - hsurfs(3,jdim-2,1)
    hsurfs(1,jdim,2) =     hsurf (1,jdim  )   - hsurfs(1,jdim  ,1)

    ! extrapolation of NE point
    hsurfs(idim,jdim,1) = 2*hsurfs(idim-1,jdim-1,1)-hsurfs(idim-2,jdim-2,1)
    hsurfs(idim,jdim,2) =   hsurf (idim  ,jdim)    -hsurfs(idim  ,jdim,1)
  ENDIF

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE sleve_split_oro

!==============================================================================
!==============================================================================
!+ Defines all subroutines for the generic routine smoother
!------------------------------------------------------------------------------
!
! SUBROUTINE smoother (finout, ie, je, nlength, nfilt)
!
!------------------------------------------------------------------------------
!
! Description:
!   This routine smoothes an arbitrary two-dimensional field (fin) by applying
!   digital filters of length nlength (4 or 8) nfilt times. The filterd field
!   is written on fout.
!
! Method:
!   Call of digital filters (dfilt4 or dfilt8) in each direction.    
!
!------------------------------------------------------------------------------
!+ Subroutine for double precision
!------------------------------------------------------------------------------

SUBROUTINE smoother_double (finout, ie, je, nlength, nfilt)

!------------------------------------------------------------------------------
!
! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)      ::    &
  ie, je,         & ! Dimension of the field
  nlength,        & ! Filter lenght
  nfilt             ! Number of iterative filerings

REAL (KIND=idouble),   INTENT (INOUT)      ::    &
  finout (ie*je)    ! 2-d field: unfiltered at input, filtered at output

! Local variables
INTEGER (KIND=iintegers) ::    &
  i,j               ! loop indicees

REAL (KIND=idouble)      ::    &
  f_2d_field(ie,je),             & !
  sxin(ie), sxh(ie), sxout(ie),  & ! local storage
  syin(je), syh(je), syout(je)     ! local storage

!------------------------------------------------------------------------------
! begin subroutine smoother_double

  f_2d_field = RESHAPE (finout, (/ie,je/))

  IF ( nlength /= 4 .AND. nlength /= 8 ) THEN
    PRINT*, ' CAUTION: Filterlength =',nlength,' not implemented'
    PRINT*, ' No filtering of output field done'
    RETURN
  ENDIF

  DO j = 1, je
    sxin(:) = f_2d_field(:,j)
    IF(nlength==4)  CALL dfilt4 ( sxin, ie, sxh, sxout, nfilt )
    IF(nlength==8)  CALL dfilt8 ( sxin, ie, sxh, sxout, nfilt )
    f_2d_field(:,j) = sxout(:)
  ENDDO
  DO i = 1, ie
    syin(:) = f_2d_field(i,:)
    IF(nlength==4)  CALL dfilt4 ( syin, je, syh, syout, nfilt )
    IF(nlength==8)  CALL dfilt8 ( syin, je, syh, syout, nfilt )
    f_2d_field(i,:) = syout(:)
  ENDDO

  finout = RESHAPE (f_2d_field, (/ie*je/))

END SUBROUTINE smoother_double

!------------------------------------------------------------------------------
!+ Subroutine for single precision
!------------------------------------------------------------------------------

SUBROUTINE smoother_single (finout, ie, je, nlength, nfilt)

!------------------------------------------------------------------------------
!
! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)      ::    &
  ie, je,         & ! Dimension of the field
  nlength,        & ! Filter lenght
  nfilt             ! Number of iterative filerings

REAL (KIND=isingle),   INTENT (INOUT)      ::    &
  finout (ie*je)    ! 2-d field: unfiltered at input, filtered at output

! Local variables
INTEGER (KIND=iintegers) ::    &
  i,j               ! loop indicees

REAL (KIND=isingle)      ::    &
  f_2d_field(ie,je),             & !
  sxin(ie), sxh(ie), sxout(ie),  & ! local storage
  syin(je), syh(je), syout(je)     ! local storage

!------------------------------------------------------------------------------
! begin subroutine smoother_single

  f_2d_field = RESHAPE (finout, (/ie,je/))

  IF ( nlength /= 4 .AND. nlength /= 8 ) THEN
    PRINT*, ' CAUTION: Filterlength =',nlength,' not implemented'
    PRINT*, ' No filtering of output field done'
    RETURN
  ENDIF

  DO j = 1, je
    sxin(:) = f_2d_field(:,j)
    IF(nlength==4)  CALL dfilt4 ( sxin, ie, sxh, sxout, nfilt )
    IF(nlength==8)  CALL dfilt8 ( sxin, ie, sxh, sxout, nfilt )
    f_2d_field(:,j) = sxout(:)
  ENDDO
  DO i = 1, ie
    syin(:) = f_2d_field(i,:)
    IF(nlength==4)  CALL dfilt4 ( syin, je, syh, syout, nfilt )
    IF(nlength==8)  CALL dfilt8 ( syin, je, syh, syout, nfilt )
    f_2d_field(i,:) = syout(:)
  ENDDO

  finout = RESHAPE (f_2d_field, (/ie*je/))

END SUBROUTINE smoother_single

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE smooth9 (finout, imin,imaxx,jmin,jmaxx,ie,je,ke)

!------------------------------------------------------------------------------
!
! Description:
!   This routine smoothes an arbitrary two-dimensional field (finout) by applying
!   a 9 points smoother. The filtered field is written on finout.
!
! Method:
!   A 9 points smoother is applied for the computation of the
!   large scale part of finout. The boundary values are treated
!   separately to assure, that also these points are smoothed.
!
!------------------------------------------------------------------------------

! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)      ::    &
  imin,jmin,imaxx,jmaxx,    & ! Local Dimension of the field
  ie, je, ke                  ! Dimension of the field

REAL (KIND=ireals), INTENT (INOUT)      ::    &
  finout (ie,je,ke)    ! 3-d field: unfiltered at input, filtered at output

! Local variables
INTEGER (KIND=iintegers) ::    &
  i,j,k,            &! loop indicees
  imin1, imaxx1,jmin1, jmaxx1

LOGICAL lbord12,lbord13,lbord24,lbord34,lcorn1,lcorn2,lcorn3,lcorn4

REAL (KIND=ireals) ::  fhelp (ie,je) ! local storage

!------------------------------------------------------------------------------
! begin subroutine smooth9
!
      lcorn1=.false.
      lcorn2=.false.
      lcorn3=.false.
      lcorn4=.false.
      lbord13=.false.
      lbord12=.false.
      lbord34=.false.
      lbord24=.false.

! If applied to a subdomain, smooth also points of a grid line halo

      imin1=imin
      jmin1=jmin
      imaxx1=imaxx
      jmaxx1=jmaxx

! Adapt dimensions to the subdomain type (if it has edge and/or corner points)

      IF( imin == 1 ) THEN
        imin1=imin+1
        lbord12=.true.
        IF(jmin == 1) THEN
          lcorn1=.true.
        ENDIF
        IF(jmaxx == je) THEN
          lcorn2=.true.
        ENDIF
      ENDIF
      IF( jmin == 1 ) THEN
        jmin1=jmin+1
        lbord13=.true.
      ENDIF
      IF(imaxx == ie )THEN
        imaxx1=imaxx-1
        lbord34=.true.
        IF (jmin == 1) THEN
          lcorn3=.true.
        ENDIF
        IF(jmaxx == je) THEN
          lcorn4=.true.
        ENDIF
      ENDIF
      IF(jmaxx == je) THEN
        jmaxx1=jmaxx-1
        lbord24=.true.
      ENDIF

      DO k = 1, ke

        fhelp(:,:) = finout(:,:,k)

! Treat inner points

        DO i=imin1,imaxx1
          DO j=jmin1,jmaxx1
            finout(i,j,k) = 0.25_ireals * fhelp(i,j)                      &
              + 0.125_ireals  * (fhelp(i-1,j  ) + fhelp(i+1,j ) + &
                                 fhelp(i  ,j-1) + fhelp(i  ,j+1))  &
              + 0.0625_ireals * (fhelp(i-1,j-1) + fhelp(i+1,j-1) + &
                                 fhelp(i-1,j+1) + fhelp(i+1,j+1))
          ENDDO
        ENDDO

! Treat corner points

        IF (lcorn1) finout(1,1,k) = 0.25_ireals *                                &
                         (fhelp(1   ,1   ) + fhelp(2   ,1   )      &
                        + fhelp(1   ,2   ) + fhelp(2   ,2   ))

        IF (lcorn2) finout(1 ,jmaxx,k) =  0.25_ireals *                           &
                         (fhelp(1   ,jmaxx  ) + fhelp(   2,jmaxx  )  &
                        + fhelp(1   ,jmaxx-1) + fhelp(   2,jmaxx-1))

        IF (lcorn3) finout(imaxx,1 ,k) =  0.25_ireals *                           &
                         (fhelp(imaxx,1  ) + fhelp(imaxx-1,1 )    &
                        + fhelp(imaxx,2  ) + fhelp(imaxx-1,2 ))

        IF (lcorn4) finout(imaxx,jmaxx,k) =  0.25_ireals *                         &
                         (fhelp(imaxx,jmaxx ) + fhelp(imaxx-1,jmaxx  )&
                        + fhelp(imaxx,jmaxx-1) + fhelp(imaxx-1,jmaxx-1))
! Treat edge points

        IF(lbord12)THEN
          DO j = jmin1,jmaxx1
            finout(1,j,k)  =                                                     &
               0.25_ireals  * (fhelp(1   ,j ) + fhelp(2   ,j ))    &
               +  0.125_ireals * (fhelp(1   ,j -1) + fhelp(1   ,j +1) +   &
                           fhelp(2   ,j -1) + fhelp(2   ,j +1) )
          ENDDO
        ENDIF

        IF(lbord34)THEN
          DO j = jmin1,jmaxx1
            finout(imaxx,j,k) =                                                   &
               0.25_ireals  * (fhelp(imaxx  ,j   ) + fhelp(imaxx-1,j  ))  &
              +  0.125_ireals * (fhelp(imaxx  ,j -1) + fhelp(imaxx  ,j +1) + &
                         fhelp(imaxx-1,j -1) + fhelp(imaxx-1,j +1) )
          ENDDO
        ENDIF

        IF(lbord13)THEN
          DO i = imin1,imaxx1
            finout(i,1,k)  =                                                     &
                  0.25_ireals  * (fhelp(i   ,1  ) + fhelp(i   ,2  ))    &
              +   0.125_ireals * (fhelp(i -1,1 ) + fhelp(i +1,1  ) +   &
                           fhelp(i -1,2 ) + fhelp(i +1,2 ) )
          ENDDO
        ENDIF

        IF(lbord24)THEN
          DO i = imin1,imaxx1
            finout(i,jmaxx,k) =                                                   &
             0.25_ireals  * (fhelp(i   ,jmaxx ) + fhelp(i   ,jmaxx-1))  &
             +  0.125_ireals * (fhelp(i -1,jmaxx ) + fhelp(i +1,jmaxx ) + &
                         fhelp(i -1,jmaxx-1) + fhelp(i +1,jmaxx-1) )
          ENDDO
        ENDIF

      ENDDO


END SUBROUTINE smooth9

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE TAUTSP(TAU,GTAU,NTAU,GAMMA,S,BREAK,COEF,L,IFLAG)

!------------------------------------------------------------------------------
!
! Description:
!   BERECHNET DEN TAUT-SPLINE FUER DIE DATEN: TAU(I),GTAU(I),I=1,.,NTAU
!
! Method:
!   WENN GAMMA GT.0  WERDEN ZUSAETZLICHE KNOTEN BERECHNET
!   GAMMA I.A.  2.5   BZW.  5.5
!
!   BREAK,COEF,L,K GEBEN DIE PP-DARSTELLUNG DER INTERPOLATIONSFUNKTION
!
!   FUER BREAK(I).LE.X.LE.BREAK(I+1) HAT DIE INTERPOLATIONSFUNKTION
!   FOLGENDE FORM
!
!   F(X)=COEF(1,I)+DX(COEF(2,I)+DX/2(COEF(3,I)+DX/3(COEF(4,I)))
!   MIT   DX=X-BREAK(I) UND I=1,...,L
!
!   IFLAG=0  KEIN FEHLER IN TAUTSP
!   IFLAG=2  FALSCHER INPUT
!
!   S(NTAU,6)  WORK-ARRAY
!
!==============================================================================

INTEGER (KIND=iintegers) IFLAG,L,NTAU,I,METHOD,NTAUM1
REAL (KIND=ireals)                                                 &
   BREAK(L),COEF(4,L),GAMMA,GTAU(NTAU),S(NTAU,6),TAU(NTAU),        &
   ALPHA,C,D,DEL,DENOM,DIVDIF,ENTRY,ENTRY3,FACTOR,FACTR2,GAM,      &
   ONEMG3,ONEMZT,RATIO,SIXTH,TEMP,X,Z,ZETA,ZT2,ALPH

!==============================================================================
!
      ALPH(X)= MIN (1.0_ireals,ONEMG3/X)
!
!   TEST DER INPUTPARAMETER
!
      IF (NTAU .LT. 4) THEN
        PRINT 600,NTAU
 600    FORMAT('  NTAU =',I4,'  NTAU SOLL GROESSER ALS 4 SEIN')
        GO TO 999
      ENDIF
!
!   BERECHNUNG VON DELTA TAU UND DER 1. UND 2.ABLEITUNG DER DATEN
!
      NTAUM1=NTAU-1
      DO I=1,NTAUM1
      S(I,1)=TAU(I+1)-TAU(I)
      IF (S(I,1) .LE. 0.0) THEN
        PRINT 610,I,TAU(I),TAU(I+1)
 610    FORMAT(' PUNKT ',I3,'  UND DIE NAECHSTEN',2E15.6,'SIND IN&
     &  FALSCHER REIHENFOLGE')
        GO TO 999
      ENDIF
      S(I+1,4) = (GTAU(I+1) - GTAU(I))/S(I,1)
      ENDDO
!
      DO I=2,NTAUM1
      S(I,4) = S(I+1,4) - S(I,4)
      ENDDO
!
!   2.ABLEITUNG VON GTAU AN DEN PUNKTEN TAU
!
      I=2
      S(2,2) = S(1,1)/3.0
      SIXTH = 1.0/6.0
      METHOD = 2
      GAM = GAMMA
      IF(GAM .LE. 0.0) METHOD = 1
      IF(GAM .GT. 3.0) THEN
        METHOD = 3
        GAM = GAM - 3.0
      ENDIF
      ONEMG3=1.0 - GAM/3.0
!
!   SCHLEIFE UEBER I
!
 10   CONTINUE
      Z=0.5
      IF (METHOD .EQ. 1) THEN
        GO TO 19
      ELSE IF (METHOD .EQ. 2) THEN
        GO TO 11
      ELSE IF (METHOD .EQ. 3) THEN
        GO TO 12
      ENDIF
 11   CONTINUE
      IF(S(I,4)*S(I+1,4).LT.0.0) GO TO 19
 12   CONTINUE
      TEMP = ABS(S(I+1,4))
      DENOM = ABS(S(I,4)) + TEMP
      IF(DENOM.EQ.0.0) GO TO 19
      Z = TEMP/DENOM
      IF(ABS(Z-0.5).LE.SIXTH) Z=0.5
 19   CONTINUE
      S(I,5) = Z
!
!   ERSTELLEN EINES TEILES DER I-TEN GLEICHUNG
!
      IF (Z-0.5 .LT. 0.) THEN
        ZETA = GAM*Z
        ONEMZT = 1.0 - ZETA
        ZT2 = ZETA**2
        ALPHA = ALPH(ONEMZT)
        FACTOR = ZETA/(ALPHA*(ZT2 - 1.0) + 1.0)
        S(I,6) = ZETA*FACTOR/6.0
        S(I,2) = S(I,2) + S(I,1)*((1.0-ALPHA*ONEMZT)*FACTOR*0.5-S(I,6))
        IF(S(I,2).LE.0.0) S(I,2) = 1.0
        S(I,3) = S(I,1)/6.0
!
      ELSE IF (Z-0.5 .EQ. 0.) THEN
!
        S(I,2) = S(I,2) + S(I,1)/3.0
        S(I,3) = S(I,1)/6.0
!
      ELSE
!
        ONEMZT = GAM*(1.0 - Z)
        ZETA = 1.0 - ONEMZT
        ALPHA = ALPH(ZETA)
        FACTOR = ONEMZT/(1.0 - ALPHA*ZETA*(1.0 + ONEMZT))
        S(I,6) = ONEMZT*FACTOR/6.0
        S(I,2) = S(I,2) + S(I,1)/3.0
        S(I,3) = S(I,6) * S(I,1)
      ENDIF
!
      IF (I .GT. 2) GO TO 30
      S(1,5) = 0.5
!
!   DIE ERSTEN BEIDEN GLEICHUNGEN ERZWINGEN STETIGKEIT DER 1. UND 3.AB-
!   LEITUNG IN TAU(I)
!
      S(1,2) = S(1,1)/6.0
      S(1,3) = S(2,2)
      ENTRY3 = S(2,3)

      IF (Z-0.5 .LT. 0.) THEN
        FACTR2 = ZETA*(ALPHA*(ZT2-1.0)+1.0)/(ALPHA*(ZETA*ZT2-1.0) + 1.0)
        RATIO = FACTR2*S(2,1)/S(1,2)
        S(2,2) = FACTR2*S(2,1) + S(1,1)
        S(2,3) = -FACTR2 * S(1,1)
!
      ELSE IF (Z-0.5 .EQ. 0.) THEN
!
        RATIO = S(2,1)/S(1,2)
        S(2,2) = S(2,1) + S(1,1)
        S(2,3) = -S(1,1)
!
      ELSE
!
        RATIO = S(2,1)/S(1,2)
        S(2,2) = S(2,1) + S(1,1)
        S(2,3) = -S(1,1)*6.0*ALPHA*S(2,6)
      ENDIF
!
!   ELIMINATION DER 1.UNBEKANNTEN AUS DER 2.GLEICHUNG
!
      S(2,2) = RATIO*S(1,3) + S(2,2)
      S(2,3) = RATIO*ENTRY3 + S(2,3)
      S(1,4) = S(2,4)
      S(2,4) = RATIO*S(1,4)
      GO TO 35
!
!
 30   CONTINUE
      S(I,2) = RATIO*S(I-1,3) + S(I,2)
      S(I,4) = RATIO*S(I-1,4) + S(I,4)
!
!   AUFSTELLEN DES TEILES DER NAECHSTEN GLEICHUNG,DER VOM I-TEN INTERVAL
!   ABHAENGT
!
 35   CONTINUE
      IF (Z-0.5 .LT. 0.) THEN
        RATIO = -S(I,6)*S(I,1)/S(I,2)
        S(I+1,2) = S(I,1)/3.0
!
      ELSE IF (Z-0.5 .EQ. 0.) THEN
!
        RATIO = -S(I,1)/(6.0*S(I,2))
        S(I+1,2) = S(I,1)/3.0
!
      ELSE
!
        RATIO = -(S(I,1)/6.0)/S(I,2)
        S(I+1,2) = S(I,1)*((1.0 - ZETA*ALPHA)*0.5*FACTOR-S(I,6))
      ENDIF
!
!   ENDE DER SCHLEIFE UEBER I
!
      I = I + 1
      IF(I.LT.NTAUM1) GO TO 10
      S(I,5) = 0.5
!
!   DIE BEIDEN LETZTEN GLEICHUNGEN ERZWINGEN STETIGKEIT DER
!   1. UND 3. ABLEITUNG IN TAU(NTAU-1)
!
      ENTRY = RATIO*S(I-1,3) + S(I,2) + S(I,1)/3.0
      S(I+1,2) = S(I,1)/6.0
      S(I+1,4) = RATIO*S(I-1,4) + S(I,4)
      IF (Z-0.5 .LT. 0.) THEN
        RATIO = S(I,1)*6.0*S(I-1,6)*ALPHA/S(I-1,2)
        S(I,2) = RATIO*S(I-1,3) + S(I,1) + S(I-1,1)
        S(I,3) = -S(I-1,1)
!
      ELSE IF (Z-0.5 .EQ. 0.) THEN
!
        RATIO = S(I,1)/S(I-1,2)
        S(I,2) = RATIO*S(I-1,3) + S(I,1) + S(I-1,1)
        S(I,3) = -S(I-1,1)
!
      ELSE
!
        FACTR2 = ONEMZT*(ALPHA*(ONEMZT**2-1.0)+1.0)/     &
                        (ALPHA*(ONEMZT**3-1.0)+1.0)
        RATIO = FACTR2*S(I,1)/S(I-1,2)
        S(I,2) = RATIO*S(I-1,3) + FACTR2*S(I-1,1) + S(I,1)
        S(I,3) = -FACTR2*S(I-1,1)
      ENDIF
!
!   ELIMINATION VON XI
!
      S(I,4) = RATIO*S(I-1,4)
      RATIO = -ENTRY/S(I,2)
      S(I+1,2) = RATIO*S(I,3) + S(I+1,2)
      S(I+1,4) = RATIO*S(I,4) + S(I+1,4)
!
!   RUECKSUBSTITUTION
!
      S(NTAU,4) = S(NTAU,4)/S(NTAU,2)
 50   CONTINUE
      S(I,4) = (S(I,4) - S(I,3)*S(I+1,4))/S(I,2)
      I = I - 1
      IF(I.GT.1) GO TO 50
      S(1,4) = (S(1,4) -S(1,3)*S(2,4)-ENTRY3*S(3,4))/S(1,2)
!
!   ERZEUGEN DER POLYNOM-TEILE
!
      BREAK(1) = TAU(1)
      L = 1
      DO 70 I = 1,NTAUM1
      COEF(1,L) = GTAU(I)
      COEF(3,L) = S(I,4)
      DIVDIF = (GTAU(I+1) - GTAU(I))/S(I,1)
      Z = S(I,5)
!
      IF (Z-0.5 .LT. 0.) THEN
        IF(Z.EQ.0.0) GO TO 65
        ZETA = GAM*Z
        ONEMZT = 1.0 -ZETA
        C = S(I+1,4)/6.0
        D = S(I,4)*S(I,6)
        L = L + 1
        DEL = ZETA*S(I,1)
        BREAK(L) = TAU(I) + DEL
        ZT2 = ZETA**2
        ALPHA = ALPH(ONEMZT)
        FACTOR = ONEMZT**2*ALPHA
        COEF(1,L) = GTAU(I) + DIVDIF*DEL+S(I,1)**2*(D*ONEMZT*(FACTOR- &
           1.0) + C*ZETA*(ZT2 - 1.0))
        COEF(2,L) = DIVDIF + S(I,1)*(D*(1.0-3.0*FACTOR) + C*(3.0*ZT2- &
           1.0))
        COEF(3,L) = 6.0*(D*ALPHA*ONEMZT + C*ZETA)
        COEF(4,L) = 6.0*(C - D*ALPHA)/S(I,1)
        COEF(4,L-1) = COEF(4,L) -6.0*D*(1.0-ALPHA)/(DEL*ZT2)
        COEF(2,L-1) = COEF(2,L) - DEL*(COEF(3,L)-DEL*0.5*COEF(4,L-1))
        GO TO 68
!
      ELSE IF (Z-0.5 .EQ. 0.) THEN
!
        COEF(2,L) = DIVDIF - S(I,1)*(2.0*S(I,4) + S(I+1,4))/6.0
        COEF(4,L) = (S(I+1,4) - S(I,4))/S(I,1)
        GO TO 68
!
      ELSE
!
        ONEMZT = GAM*(1.0 - Z)
        IF(ONEMZT.EQ.0.0) GO TO 65
        ZETA = 1.0 - ONEMZT
        ALPHA = ALPH(ZETA)
        C = S(I+1,4)*S(I,6)
        D = S(I,4)/6.0
        DEL = ZETA*S(I,1)
        BREAK(L+1) = TAU(I) + DEL
        COEF(2,L) = DIVDIF -S(I,1)*(2.0*D + C)
        COEF(4,L) = 6.0*(C*ALPHA - D)/S(I,1)
        L = L + 1
        COEF(4,L) = COEF(4,L-1) + 6.0*(1.0-ALPHA)*C/(S(I,1)*ONEMZT**3)
        COEF(3,L) = COEF(3,L-1) + DEL* COEF(4,L-1)
        COEF(2,L) = COEF(2,L-1) + DEL*(COEF(3,L-1)+DEL*0.5*COEF(4,L-1))
        COEF(1,L) = COEF(1,L-1) + DEL*(COEF(2,L-1)+DEL*0.5*   &
              (COEF(3,L-1) + (DEL/3.0)*COEF(4,L-1)))
        GO TO 68
      ENDIF
!
!
   65 CONTINUE
      COEF(2,L) = DIVDIF
      COEF(3,L) = 0.0
      COEF(4,L) = 0.0
   68 CONTINUE
!
      L = L + 1
      BREAK(L) = TAU(I+1)
   70 CONTINUE

      IFLAG = 0
      RETURN
!
 999  IFLAG = 2

END SUBROUTINE tautsp

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE tautsp2D (TAU, GTAU, NTAU, NI, IMIN, IMAX, NTAUMAX, GAMMA,       &
                     S, BREAK, COEF, L_VEC, IFLAG)

!------------------------------------------------------------------------------
!
! Description:
!   BERECHNET DEN TAUT-SPLINE FUER DIE DATEN: TAU(I),GTAU(I),I=1,.,NTAU
!
! Method:
!   WENN GAMMA GT.0  WERDEN ZUSAETZLICHE KNOTEN BERECHNET
!   GAMMA I.A.  2.5   BZW.  5.5
!
!   BREAK,COEF,L,K GEBEN DIE PP-DARSTELLUNG DER INTERPOLATIONSFUNKTION
!
!   FUER BREAK(I).LE.X.LE.BREAK(I+1) HAT DIE INTERPOLATIONSFUNKTION
!   FOLGENDE FORM
!
!   F(X)=COEF(1,I)+DX(COEF(2,I)+DX/2(COEF(3,I)+DX/3(COEF(4,I)))
!   MIT   DX=X-BREAK(I) UND I=1,...,L
!
!   IFLAG=0  KEIN FEHLER IN TAUTSP
!   IFLAG=2  FALSCHER INPUT
!
!   S(NTAU,6)  WORK-ARRAY
!
!==============================================================================

INTEGER(KIND=iintegers), INTENT(IN)   :: &
    NI, IMIN, IMAX, NTAUMAX

INTEGER(KIND=iintegers), INTENT(IN)   :: &
    NTAU(NI)

REAL   (KIND=ireals),    INTENT(IN)   :: &
    GTAU(NI,NTAUMAX),                    &
    TAU (NI,NTAUMAX),                    &
    GAMMA

REAL   (KIND=ireals),    INTENT(OUT)  :: &
    BREAK(NI,*),                         &
    COEF (NI,4,*),                       &
    S    (NI,NTAUMAX,6)

INTEGER(KIND=iintegers), INTENT(OUT)  :: &
    L_VEC(NI),                           &
    IFLAG

! Local Variables
INTEGER(KIND=iintegers)  I,J,K,L,METHOD,NTAUM1, mb_err_indx_i,mb_err_indx_k

REAL(KIND=ireals) C,D,DEL,DENOM,DIVDIF,ENTRY,FACTR2,GAM,      &
                  ONEMG3,SIXTH,TEMP,X,ALPH

REAL(KIND=ireals)                     :: &
    RATIO_VEC   (NI),                    &
    Z_VEC       (NI),                    &
    ZETA_VEC    (NI),                    &
    ZT2_VEC     (NI),                    &
    ALPHA_VEC   (NI),                    &
    FACTOR_VEC  (NI),                    &
    ONEMZT_VEC  (NI),                    &
    ENTRY3      (NI)

!==============================================================================
!
   ALPH(X)= MIN (1.0_ireals,ONEMG3/X)
!
!   TEST DER INPUTPARAMETER
!
      mb_err_indx_i = -1
      DO i = IMIN, IMAX
         IF (NTAU(i) .LT. 4) mb_err_indx_i = i
      ENDDO

      IF ( mb_err_indx_i /= -1 ) THEN
        PRINT *, '  NTAU =', NTAU(mb_err_indx_i), mb_err_indx_i, &
                 ':  MUST BE BIGGER THAN 4'
!        PRINT 600,NTAU(i)
!600     FORMAT('  NTAU =',I4,'  NTAU SOLL GROESSER ALS 4 SEIN')
         GO TO 999
      ENDIF

!
!   BERECHNUNG VON DELTA TAU UND DER 1. UND 2.ABLEITUNG DER DATEN
!
      mb_err_indx_i = -1
      mb_err_indx_k = -1

      DO k = 1, NTAUMAX
         DO i = IMIN, IMAX
            IF (k <= NTAU(i)-1) THEN
               S(i,k,1)=TAU(i,k+1)-TAU(i,k)
               IF (S(i,k,1) .LE. 0.0) THEN
                  mb_err_indx_i = i
                  mb_err_indx_k = k
               ELSE
                  S(i,k+1,4) = (GTAU(i,k+1) - GTAU(i,k))/S(i,k,1)
               ENDIF
            ENDIF
         ENDDO

         IF ( mb_err_indx_i /= -1 .OR. mb_err_indx_k /= -1) THEN
            PRINT 610,mb_err_indx_i,mb_err_indx_k, &
     &                TAU(mb_err_indx_i,mb_err_indx_k),&
     &                TAU(mb_err_indx_i,mb_err_indx_k+1)
 610        FORMAT(' PUNKT ',2I3,'  UND DIE NAECHSTEN',2E15.6,'SIND IN&
     &      FALSCHER REIHENFOLGE')
            GO TO 999
         ENDIF
      ENDDO

      DO k = 1, NTAUMAX
         DO i = IMIN, IMAX
            IF (k >= 2 .AND. k <= NTAU(i)-1) THEN
               S(i,k,4) = S(i,k+1,4) - S(i,k,4)
            ENDIF
         ENDDO
      ENDDO
!
!   2.ABLEITUNG VON GTAU AN DEN PUNKTEN TAU
!
      DO i = IMIN, IMAX
         S(i,2,2) = S(i,1,1)/3.0
      ENDDO

      SIXTH = 1.0/6.0
      METHOD = 2
      GAM = GAMMA
      IF(GAM .LE. 0.0) METHOD = 1
      IF(GAM .GT. 3.0) THEN
        METHOD = 3
        GAM = GAM - 3.0
      ENDIF
      ONEMG3=1.0 - GAM/3.0
!
!   SCHLEIFE UEBER K
!
      DO k = 2, NTAUMAX
         DO i = IMIN, IMAX
            IF ( k <= NTAU(i)-2 ) THEN

               Z_VEC(i)=0.5
               IF (METHOD /= 1) THEN
                  IF ( ((METHOD == 2) .AND.                                   &
                    (S(i,k,4)*S(i,k+1,4) >= 0.0)) .OR. (METHOD == 3) ) THEN
                     TEMP = ABS(S(i,k+1,4))
                     DENOM = ABS(S(i,k,4)) + TEMP
                     IF (DENOM /= 0.0) THEN
                        Z_VEC(i) = TEMP/DENOM
                        IF(ABS(Z_VEC(i)-0.5).LE.SIXTH) Z_VEC(i)=0.5
                     ENDIF
                  ENDIF
               ENDIF
               S(i,k,5) = Z_VEC(i)
!
!   ERSTELLEN EINES TEILES DER I-TEN GLEICHUNG
!
               IF (Z_VEC(i)-0.5 .LT. 0.) THEN
                  ZETA_VEC(i) = GAM*Z_VEC(i)
                  ONEMZT_VEC(i) = 1.0 - ZETA_VEC(i)
                  ZT2_VEC(i) = ZETA_VEC(i)**2
                  ALPHA_VEC(i) = ALPH(ONEMZT_VEC(i))
                  FACTOR_VEC(i) = ZETA_VEC(i) /                       &
                                   (ALPHA_VEC(i)*(ZT2_VEC(i) - 1.0) + 1.0)
                  S(i,k,6) = ZETA_VEC(i)*FACTOR_VEC(i)/6.0
                  S(i,k,2) = S(i,k,2) + S(i,k,1) *                  &
                                   ((1.0-ALPHA_VEC(i)*ONEMZT_VEC(i))  &
                                     *FACTOR_VEC(i)*0.5-S(i,k,6))
                  IF(S(i,k,2).LE.0.0) S(i,k,2) = 1.0
                  S(i,k,3) = S(i,k,1)/6.0
!
               ELSE IF (Z_VEC(i)-0.5 .EQ. 0.) THEN
!
                  S(i,k,2) = S(i,k,2) + S(i,k,1)/3.0
                  S(i,k,3) = S(i,k,1)/6.0
!
               ELSE
!
                  ONEMZT_VEC(i) = GAM*(1.0 - Z_VEC(i))
                  ZETA_VEC(i) = 1.0 - ONEMZT_VEC(i)
                  ALPHA_VEC(i) = ALPH(ZETA_VEC(i))
                  FACTOR_VEC(i) = ONEMZT_VEC(i) /                    &
                   (1.0 - ALPHA_VEC(i)*ZETA_VEC(i)*(1.0 + ONEMZT_VEC(i)))
                  S(i,k,6) = ONEMZT_VEC(i)*FACTOR_VEC(i)/6.0
                  S(i,k,2) = S(i,k,2) + S(i,k,1)/3.0
                  S(i,k,3) = S(i,k,6) * S(i,k,1)
               ENDIF
!
               IF (k == 2) THEN
                  S(i,1,5) = 0.5
!
!   DIE ERSTEN BEIDEN GLEICHUNGEN ERZWINGEN STETIGKEIT DER 1. UND 3.AB-
!   LEITUNG IN TAU(i,k)
!
                  S(i,1,2) = S(i,1,1)/6.0
                  S(i,1,3) = S(i,2,2)
                  ENTRY3(i) = S(i,2,3)

                  IF (Z_VEC(i)-0.5 .LT. 0.) THEN
                     FACTR2 = ZETA_VEC(i)*(ALPHA_VEC(i)                 &
                              *(ZT2_VEC(i)-1.0)+1.0)/                     &
                        (ALPHA_VEC(i)*(ZETA_VEC(i)*ZT2_VEC(i)-1.0) + 1.0)
                     RATIO_VEC(i) = FACTR2*S(i,2,1)/S(i,1,2)
                     S(i,2,2) = FACTR2*S(i,2,1) + S(i,1,1)
                     S(i,2,3) = -FACTR2 * S(i,1,1)
!
                  ELSE IF (Z_VEC(i)-0.5 .EQ. 0.) THEN
!
                     RATIO_VEC(i) = S(i,2,1)/S(i,1,2)
                     S(i,2,2) = S(i,2,1) + S(i,1,1)
                     S(i,2,3) = -S(i,1,1)
!
                  ELSE
!
                     RATIO_VEC(i) = S(i,2,1)/S(i,1,2)
                     S(i,2,2) = S(i,2,1) + S(i,1,1)
                     S(i,2,3) = -S(i,1,1)*6.0*ALPHA_VEC(i)*S(i,2,6)
                  ENDIF
!
!   ELIMINATION DER 1.UNBEKANNTEN AUS DER 2.GLEICHUNG
!
                  S(i,2,2) = RATIO_VEC(i)*S(i,1,3) + S(i,2,2)
                  S(i,2,3) = RATIO_VEC(i)*ENTRY3(i) + S(i,2,3)
                  S(i,1,4) = S(i,2,4)
                  S(i,2,4) = RATIO_VEC(i)*S(i,1,4)
!

               ELSE ! k > 2
!
                  S(i,k,2) = RATIO_VEC(i)*S(i,k-1,3) + S(i,k,2)
                  S(i,k,4) = RATIO_VEC(i)*S(i,k-1,4) + S(i,k,4)
               ENDIF  ! k == 2
!
!   AUFSTELLEN DES TEILES DER NAECHSTEN GLEICHUNG,DER VOM I-TEN INTERVAL
!   ABHAENGT
!
               IF (Z_VEC(i)-0.5 .LT. 0.) THEN
                  RATIO_VEC(i) = -S(i,k,6)*S(i,k,1)/S(i,k,2)
                  S(i,k+1,2) = S(i,k,1)/3.0
!
               ELSE IF (Z_VEC(i)-0.5 .EQ. 0.) THEN
!
                  RATIO_VEC(i) = -S(i,k,1)/(6.0*S(i,k,2))
                  S(i,k+1,2) = S(i,k,1)/3.0
!
               ELSE
!
                  RATIO_VEC(i) = -(S(i,k,1)/6.0)/S(i,k,2)
                  S(i,k+1,2)   = S(i,k,1) *                            &
                               ((1.0 - ZETA_VEC(i)*ALPHA_VEC(i))*      &
                               0.5*FACTOR_VEC(i)-S(i,k,6))
               ENDIF
!
!   ENDE DER SCHLEIFE UEBER k
!
            ENDIF ! k <= NTAU(i)-2
         ENDDO ! i = IMIN, IMAX
      ENDDO ! k = 2, NTAUMAX

      DO i = IMIN, IMAX
         k = NTAU(i)-1
         S(i,k,5) = 0.5

!
!   DIE BEIDEN LETZTEN GLEICHUNGEN ERZWINGEN STETIGKEIT DER
!   1. UND 3. ABLEITUNG IN TAU(NTAU-1)
!
         ENTRY = RATIO_VEC(i)*S(i,k-1,3) + S(i,k,2) + S(i,k,1)/3.0
         S(i,k+1,2) = S(i,k,1)/6.0
         S(i,k+1,4) = RATIO_VEC(i)*S(i,k-1,4) + S(i,k,4)
         IF (Z_VEC(i)-0.5 .LT. 0.) THEN
            RATIO_VEC(i) = S(i,k,1) * 6.0 * S(i,k-1,6) *              &
                             ALPHA_VEC(i)/S(i,k-1,2)
            S(i,k,2) = RATIO_VEC(i)*S(i,k-1,3) +S(i,k,1) + S(i,k-1,1)
            S(i,k,3) = -S(i,k-1,1)
!
         ELSE IF (Z_VEC(i)-0.5 .EQ. 0.) THEN
!
            RATIO_VEC(i) = S(i,k,1)/S(i,k-1,2)
            S(i,k,2) = RATIO_VEC(i)*S(i,k-1,3) + S(i,k,1) + S(i,k-1,1)
            S(i,k,3) = -S(i,k-1,1)
!
         ELSE
!
            FACTR2 = ONEMZT_VEC(i) * (ALPHA_VEC(i) *                    &
                            (ONEMZT_VEC(i)**2-1.0)+1.0)  /                &
                            (ALPHA_VEC(i)*(ONEMZT_VEC(i)**3-1.0)+1.0)
            RATIO_VEC(i) = FACTR2*S(i,k,1)/S(i,k-1,2)
            S(i,k,2) = RATIO_VEC(i)*S(i,k-1,3) + FACTR2*S(i,k-1,1)  &
                         + S(i,k,1)
            S(i,k,3) = -FACTR2*S(i,k-1,1)
         ENDIF
!
!   ELIMINATION VON XI
!
         S(i,k,4) = RATIO_VEC(i)*S(i,k-1,4)
         RATIO_VEC(i) = -ENTRY/S(i,k,2)
         S(i,k+1,2) = RATIO_VEC(i)*S(i,k,3) + S(i,k+1,2)
         S(i,k+1,4) = RATIO_VEC(i)*S(i,k,4) + S(i,k+1,4)
      ENDDO ! i = IMIN, IMAX

!
!   RUECKSUBSTITUTION
!
      DO i = IMIN, IMAX
         S(i,NTAU(i),4) = S(i,NTAU(i),4)/S(i,NTAU(i),2)
      ENDDO

      DO k = NTAUMAX,2,-1
         DO i = IMIN, IMAX
            IF (k <= NTAU(i)-1) THEN
               S(i,k,4) = (S(i,k,4) - S(i,k,3)*S(i,k+1,4))/S(i,k,2)
            ENDIF
         ENDDO
      ENDDO

      DO i = IMIN, IMAX
         S(i,1,4) = (S(i,1,4) - S(i,1,3)*S(i,2,4) -              &
                       ENTRY3(i)*S(i,3,4))/S(i,1,2)
      ENDDO
!
!   ERZEUGEN DER POLYNOM-TEILE
!
      DO i = IMIN, IMAX
         BREAK(i,1) = TAU(i,1)
         L_VEC(i) = 1
      ENDDO

      DO k = 1, NTAUMAX
         DO i = IMIN, IMAX
            IF ( k <= NTAU(i)-1) THEN
               L = L_VEC(i)
               COEF(i,1,L) = GTAU(i,k)
               COEF(i,3,L) = S(i,k,4)
               DIVDIF = (GTAU(i,k+1) - GTAU(i,k))/S(i,k,1)
               Z_VEC(i) = S(i,k,5)
!
               IF (Z_VEC(i)-0.5 .LT. 0.) THEN
                  ! US avoid division by 0, if Z_VEC is veeeery small
                  ! by treating it as 0
                  IF(ABS(Z_VEC(i)) < 1E-50_ireals) THEN
                     COEF(i,2,L) = DIVDIF
                     COEF(i,3,L) = 0.0
                     COEF(i,4,L) = 0.0
                  ELSE
                     ZETA_VEC(i) = GAM*Z_VEC(i)
                     ONEMZT_VEC(i) = 1.0 -ZETA_VEC(i)
                     C = S(i,k+1,4)/6.0
                     D = S(i,k,4)*S(i,k,6)
                     L = L + 1
                     DEL = ZETA_VEC(i)*S(i,k,1)
                     BREAK(i,L) = TAU(i,k) + DEL
                     ZT2_VEC(i) = ZETA_VEC(i)**2
!                    ALPHA_VEC(i) = ALPH(ONEMZT_VEC(i))
                     ALPHA_VEC(i) = MIN(1.0_ireals,ONEMG3/ONEMZT_VEC(i))
                     FACTOR_VEC(i) = ONEMZT_VEC(i)**2*ALPHA_VEC(i)
                     COEF(i,1,L) = GTAU(i,k) + DIVDIF*DEL +              &
                          S(i,k,1)**2*(D*ONEMZT_VEC(i)*(FACTOR_VEC(i)- &
                             1.0) + C*ZETA_VEC(i)*(ZT2_VEC(i) - 1.0))
                     COEF(i,2,L) = DIVDIF + S(i,k,1) *                   &
                         (D*(1.0-3.0*FACTOR_VEC(i)) + C*(3.0*ZT2_VEC(i)- &
                                    1.0))
                     COEF(i,3,L) = 6.0*(D*ALPHA_VEC(i)*ONEMZT_VEC(i)   &
                                      + C*ZETA_VEC(i))
                     COEF(i,4,L) = 6.0*(C - D*ALPHA_VEC(i))/S(i,k,1)
                     COEF(i,4,L-1) = COEF(i,4,L) -6.0*D*                 &
                                       (1.0-ALPHA_VEC(i))/(DEL*ZT2_VEC(i))
                     COEF(i,2,L-1) = COEF(i,2,L) -                       &
                                   DEL*(COEF(i,3,L)-DEL*0.5*COEF(i,4,L-1))
                  ENDIF
!
               ELSE IF (Z_VEC(i)-0.5 .EQ. 0.) THEN
!
                  COEF(i,2,L) = DIVDIF - S(i,k,1) *                      &
                                  (2.0*S(i,k,4) + S(i,k+1,4))/6.0
                  COEF(i,4,L) = (S(i,k+1,4) - S(i,k,4))/S(i,k,1)
!
               ELSE
!
                  ONEMZT_VEC(i) = GAM*(1.0 - Z_VEC(i))
                  IF(ONEMZT_VEC(i).EQ.0.0) THEN
                     COEF(i,2,L) = DIVDIF
                     COEF(i,3,L) = 0.0
                     COEF(i,4,L) = 0.0
                  ELSE
                     ZETA_VEC(i) = 1.0 - ONEMZT_VEC(i)
!                    ALPHA_VEC(i) = ALPH(ZETA_VEC(i))
                     ALPHA_VEC(i) = MIN(1.0_ireals,ONEMG3/ZETA_VEC(i))
                     C = S(i,k+1,4)*S(i,k,6)
                     D = S(i,k,4)/6.0
                     DEL = ZETA_VEC(i)*S(i,k,1)
                     BREAK(i,L+1) = TAU(i,k) + DEL
                     COEF(i,2,L) = DIVDIF -S(i,k,1)*(2.0*D + C)
                     COEF(i,4,L) = 6.0*(C*ALPHA_VEC(i) - D)/S(i,k,1)
                     L = L + 1
                     COEF(i,4,L) = COEF(i,4,L-1) + 6.0 *                 &
                         (1.0-ALPHA_VEC(i))*C/(S(i,k,1)*ONEMZT_VEC(i)**3)
                     COEF(i,3,L) = COEF(i,3,L-1) + DEL* COEF(i,4,L-1)
                     COEF(i,2,L) = COEF(i,2,L-1) + DEL*(COEF(i,3,L-1)  &
                                     +DEL*0.5*COEF(i,4,L-1))
                     COEF(i,1,L) = COEF(i,1,L-1) + DEL*(COEF(i,2,L-1)  &
                                     +DEL*0.5*                               &
                           (COEF(i,3,L-1) + (DEL/3.0)*COEF(i,4,L-1)))
                  ENDIF
               ENDIF
!
               L = L + 1
               BREAK(i,L) = TAU(i,k+1)
               L_VEC(i) = L
            ENDIF
         ENDDO ! i = IMIN, IMAX
      ENDDO ! k = 1, NTAUMAX

      IFLAG = 0

      RETURN
!
 999  IFLAG = 2

END SUBROUTINE tautsp2D

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE to_upper ( string )

!-------------------------------------------------------------------------------
!
! Description:
!   Convert alphabetic characters in 'string' from lower to upper case
!-------------------------------------------------------------------------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  CHARACTER (LEN=*), INTENT(INOUT) :: string

! Local parameters:
! ----------------
  CHARACTER (LEN=26), PARAMETER :: UPPER="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  CHARACTER (LEN=26), PARAMETER :: lower="abcdefghijklmnopqrstuvwxyz"

! Local variables:
! ---------------
  INTEGER (KIND=iintegers)      :: i, j
!
!------------ End of header ----------------------------------------------------

  DO i = 1, LEN_TRIM(string)
    j = INDEX ( lower, string(i:i) )
    IF ( j > 0 ) string(i:i) = UPPER(j:j)
  END DO

END SUBROUTINE to_upper


!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE uvrot2uv (urot, vrot, rlat, rlon, pollat, pollon, u, v)


!------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the rotated system
!   to the real geographical system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

! Parameter list:
REAL (KIND=ireals), INTENT (IN)          ::    &
  urot, vrot,     & ! wind components in the rotated grid
  rlat, rlon,     & ! latitude and longitude in the true geographical system
  pollat, pollon    ! latitude and longitude of the north pole of the
                    ! rotated grid

REAL (KIND=ireals), INTENT (OUT)         ::    &
  u, v              ! wind components in the true geographical system

! Local variables

REAL (KIND=ireals)                       ::    &
  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm

REAL (KIND=ireals)                       ::    &
  zrpi18 = 57.2957795_ireals,       & !
  zpir18 = 0.0174532925_ireals

!------------------------------------------------------------------------------
! Begin subroutine uvrot2uv
!------------------------------------------------------------------------------

! Converting from degree to radians
  zsinpol = SIN(pollat * zpir18)
  zcospol = COS(pollat * zpir18)
  zlonp   = (pollon-rlon) * zpir18
  zlat    =         rlat  * zpir18

  zarg1   = zcospol*SIN(zlonp)
  zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
  znorm   = 1.0/SQRT(zarg1**2 + zarg2**2)

! Convert the u- and v-components
  u       =   urot*zarg2*znorm + vrot*zarg1*znorm
  v       = - urot*zarg1*znorm + vrot*zarg2*znorm

END SUBROUTINE uvrot2uv

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE uvrot2uv_vec(u, v, rlat, rlon, pollat, pollon, idim, jdim)

!------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the rotated
!   system to the real geographical system. This is the vectorized form
!   of the routine above, i.e. the computation is for a whole 2D field.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Parameter list:
INTEGER (KIND=iintegers), INTENT(IN)     ::    &
  idim, jdim        ! dimensions of the field

REAL (KIND=ireals), INTENT (INOUT)       ::    &
  u  (idim,jdim), & ! wind components in the true geographical system
  v  (idim,jdim)    !

REAL (KIND=ireals), INTENT (IN)          ::    &
  rlat(idim,jdim),& ! coordinates in the true geographical system
  rlon(idim,jdim),& !
  pollat, pollon    ! latitude and longitude of the north pole of the
                    ! rotated grid

! Local variables
REAL (KIND=ireals)                       ::    &
  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm, zugeo, zvgeo

INTEGER (KIND=iintegers)                 ::    i, j
REAL (KIND=ireals)                       ::    &
  zrpi18 = 57.2957795_ireals,       & !
  zpir18 = 0.0174532925_ireals

!------------------------------------------------------------------------------
! Begin Subroutine uvrot2uv_vec
!------------------------------------------------------------------------------

! Converting from degree to radians
  zsinpol = SIN(pollat * zpir18)
  zcospol = COS(pollat * zpir18)

  DO j = 1, jdim
    DO i = 1, idim

      zlonp   = (pollon-rlon(i,j)) * zpir18
      zlat    =         rlat(i,j)  * zpir18

      zarg1   = zcospol*SIN(zlonp)
      zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
      znorm   = 1.0/SQRT(zarg1**2 + zarg2**2)

      ! Convert the u- and v-components
      zugeo   =  u(i,j)*zarg2*znorm + v(i,j)*zarg1*znorm
      zvgeo   = -u(i,j)*zarg1*znorm + v(i,j)*zarg2*znorm
      u(i,j) = zugeo
      v(i,j) = zvgeo

    ENDDO
  ENDDO

END SUBROUTINE uvrot2uv_vec

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE uv2uvrot(u, v, rlat, rlon, pollat, pollon, urot, vrot)

!------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the real
!   geographical system to the rotated system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Parameter list:
REAL (KIND=ireals), INTENT (IN)          ::    &
  u   , v   ,     & ! wind components in the true geographical system
  rlat, rlon,     & ! coordinates in the true geographical system
  pollat, pollon    ! latitude and longitude of the north pole of the
                    ! rotated grid

REAL (KIND=ireals), INTENT (OUT)         ::    &
  urot, vrot        ! wind components in the rotated grid             

! Local variables

REAL (KIND=ireals)                       ::    &
  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm

REAL (KIND=ireals)                       ::    &
  zrpi18 = 57.2957795_ireals,       & !
  zpir18 = 0.0174532925_ireals

!------------------------------------------------------------------------------
! Begin Subroutine uv2uvrot
!------------------------------------------------------------------------------

  zsinpol = SIN(pollat * zpir18)
  zcospol = COS(pollat * zpir18)
  zlonp   = (pollon-rlon) * zpir18
  zlat    =         rlat  * zpir18

  zarg1   = zcospol*SIN(zlonp)
  zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
  znorm   = 1.0_ireals/SQRT( zarg1**2 + zarg2**2 )

! Transform the u and v wind components
  urot   =  u*zarg2*znorm - v*zarg1*znorm
  vrot   =  u*zarg1*znorm + v*zarg2*znorm

END SUBROUTINE uv2uvrot

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE uv2uvrot_vec(u, v, rlat, rlon, pollat, pollon, idim, jdim)

!------------------------------------------------------------------------------
!
! Description:
!   This routine converts the wind components u and v from the real
!   geographical system to the rotated system. This is the vectorized form
!   of the routine above, i.e. the computation is for a whole 2D field.
!
! Method:
!   Transformation formulas for converting between these two systems.
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Parameter list:
INTEGER (KIND=iintegers), INTENT(IN)     ::    &
  idim, jdim        ! dimensions of the field

REAL (KIND=ireals), INTENT (INOUT)       ::    &
  u  (idim,jdim), & ! wind components in the true geographical system
  v  (idim,jdim)    !

REAL (KIND=ireals), INTENT (IN)          ::    &
  rlat(idim,jdim),& ! coordinates in the true geographical system
  rlon(idim,jdim),& !
  pollat, pollon    ! latitude and longitude of the north pole of the
                    ! rotated grid

! Local variables
REAL (KIND=ireals)                       ::    &
  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm, zurot, zvrot

INTEGER (KIND=iintegers)                 ::    i, j
REAL (KIND=ireals)                       ::    &
  zrpi18 = 57.2957795_ireals,       & !
  zpir18 = 0.0174532925_ireals

!------------------------------------------------------------------------------
! Begin Subroutine uv2uvrot_vec
!------------------------------------------------------------------------------

  zsinpol = SIN ( pollat * zpir18 )
  zcospol = COS ( pollat * zpir18 )

  DO j = 1, jdim
    DO i = 1, idim

      zlonp   = ( pollon - rlon(i,j) ) * zpir18
      zlat    =            rlat(i,j)   * zpir18

      zarg1   = zcospol*SIN(zlonp)
      zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
      znorm   = 1.0_ireals/SQRT( zarg1**2 + zarg2**2 )

      ! Transform the u and v wind components
      zurot =  u(i,j)*zarg2*znorm - v(i,j)*zarg1*znorm
      zvrot =  u(i,j)*zarg1*znorm + v(i,j)*zarg2*znorm
      u(i,j) = zurot
      v(i,j) = zvrot

    ENDDO
  ENDDO

END SUBROUTINE uv2uvrot_vec

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE uv2df (u, v, d, f)

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes wind speed and wind direction from the wind
!   components.
!
! Method:
!   Straightforward.
!
!------------------------------------------------------------------------------
!
! Parameter list:
REAL (KIND=ireals), INTENT (IN)          ::    &
  u   , v           ! wind components in the true geographical system

REAL (KIND=ireals), INTENT (OUT)         ::    &
  f   , d           ! wind speed and wind direction

! Local variables

REAL (KIND=ireals)                       ::    &
  zrpi18 = 57.2957795_ireals,       & ! conversion from radians to degrees
  zsmall = 0.001_ireals

!------------------------------------------------------------------------------
! Begin Subroutine uv2df
!------------------------------------------------------------------------------

  IF (ABS(u) > zsmall) THEN
    f  =  SQRT( u*u + v*v )
    d  =  v / u
    d  =  180.0_ireals + SIGN( 90.0_ireals , u ) - ATAN( d ) *zrpi18
  ELSEIF (ABS(v) > zsmall) THEN
    f  =  ABS( v )
    d  =  270.0_ireals - SIGN( 90.0_ireals , v )
  ELSE
    f  =    0.0_ireals
    d  =    0.0_ireals
  ENDIF

END SUBROUTINE uv2df

!==============================================================================
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE uv2df_vec (u, v, d, f, idim, jdim)

!------------------------------------------------------------------------------
!
! Description:
!   This routine computes wind speed and wind direction from the wind
!   components.  This is the vectorized form of the routine above,
!   i.e. the computation is for a whole 2D field.
!
! Method:
!   Straightforward.
!
!------------------------------------------------------------------------------
!
! Parameter list:
INTEGER (KIND=iintegers), INTENT(IN)     ::    &
  idim, jdim        ! dimensions of the field

REAL (KIND=ireals), INTENT (IN)          ::    &
  u  (idim,jdim) ,& ! wind components in the true geographical system
  v  (idim,jdim)    !

REAL (KIND=ireals), INTENT (OUT)         ::    &
  f  (idim,jdim) ,& ! wind speed
  d  (idim,jdim)    ! wind direction

! Local variables

INTEGER (KIND=iintegers)                 ::    i, j
REAL (KIND=ireals)                       ::    &
  zrpi18 = 57.2957795_ireals,       & ! conversion from radians to degrees
  zsmall = 0.001_ireals

!------------------------------------------------------------------------------
! Begin Subroutine uv2df_vec
!------------------------------------------------------------------------------

  DO j = 1, jdim
    DO i = 1, idim

      IF (ABS(u(i,j)) > zsmall) THEN
        f (i,j)  =  SQRT( u(i,j)*u(i,j) + v(i,j)*v(i,j) )
        d (i,j)  =  180.0_ireals + SIGN( 90.0_ireals , u(i,j) )               &
                                 - ATAN( v(i,j) / u(i,j) ) *zrpi18
      ELSEIF (ABS(v(i,j)) > zsmall) THEN
        f (i,j)  =  ABS( v(i,j) )
        d (i,j)  =  270.0_ireals - SIGN( 90.0_ireals , v(i,j) )
      ELSE
        f (i,j)  =    0.0_ireals
        d (i,j)  =    0.0_ireals
      ENDIF

    ENDDO
  ENDDO

END SUBROUTINE uv2df_vec

!==============================================================================
!==============================================================================

END MODULE utilities
