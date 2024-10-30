!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SFX_LIC for details. version 1.
!
MODULE MODE_SSGF
!
!
INTERFACE F09A
  MODULE PROCEDURE F09A
END INTERFACE
INTERFACE F09B
  MODULE PROCEDURE F09B
END INTERFACE
INTERFACE F09C
  MODULE PROCEDURE F09C
END INTERFACE
INTERFACE F94A
  MODULE PROCEDURE F94A
END INTERFACE
INTERFACE S20A
  MODULE PROCEDURE S20A
END INTERFACE
INTERFACE S20B
  MODULE PROCEDURE S20B
END INTERFACE
INTERFACE B22A
  MODULE PROCEDURE B22A
END INTERFACE
INTERFACE B22B
  MODULE PROCEDURE B22B
END INTERFACE
INTERFACE B21A
  MODULE PROCEDURE B21A
END INTERFACE
INTERFACE B21B
  MODULE PROCEDURE B21B
END INTERFACE
INTERFACE T18A
  MODULE PROCEDURE T18A
END INTERFACE
!
!
CONTAINS
!
!
!---------------------------------------------------------------------------------------
!
!#######################################################################################
SUBROUTINE F09A(PPHIOC,PUSR,PVMOD,PUREF,PZOSEA,PHS,PCP,PSV,PSA,PVFM,PFM)
!#######################################################################################
!
!****  *F09A*
!
!       PURPOSE
!       -------
!       Computes key parameters of the sea-spray generation SUBROUTINE.
!       BASED ON BAO's FORTRAN CODE SHARED BY FAIRALL
!       NB: PARAMS OF INTEGRATED SSGF STATS NOT GOING THROUGH ANY SIZE DISTRIBUTION
!         
!       Input Arguments:
!
!       Return:
!               PSV      --  normalized source volume (m/s)
!               PSA     --  normalized source area (1/s)
!               PVFM    --  estimate mean fall velocity (m/s)
!               PFM     --  droplet mass flux (kg/m2/s)
!
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       Fairall, C. W., Banner, M. L., Peirson, W. L., Asher, W., and Morison, R. P., 2009: 
!          Investigation of the physical scaling of sea spray spume droplet production, 
!          Journal of Geophysical Research, 114, C10 001, https://doi.org/10.1029/2008JC004918
!
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!       MODIFICATIONS
!       -------------
!---------------------------------------------------------------------------------------
!
USE MODD_CSTS,       ONLY : XG
USE MODD_OCEAN_CSTS, ONLY: XRHOSW
!
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
IMPLICIT NONE
!
!* 0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN) :: PPHIOC     ! 
REAL, DIMENSION(:), INTENT(IN) :: PUSR
REAL, DIMENSION(:), INTENT(IN) :: PVMOD
REAL, DIMENSION(:), INTENT(IN) :: PUREF
REAL, DIMENSION(:), INTENT(IN) :: PZOSEA
REAL, DIMENSION(:), INTENT(IN) :: PHS               ! Significant wave height (m)
REAL, DIMENSION(:), INTENT(IN) :: PCP               ! wave peak phase velocity (m/s)
!
REAL, DIMENSION(:), INTENT(OUT) :: PFM    ! Spray mass flux (kg m-2 s-1)
REAL, DIMENSION(:), INTENT(OUT) :: PSV    ! normalized source volume m/s NB not used subsequently for now
REAL, DIMENSION(:), INTENT(OUT) :: PSA    ! normalized source area 1/s
REAL, DIMENSION(:), INTENT(OUT) :: PVFM   ! estimate mean fall velocity (m/s)
!
!* 0.2 declarations of local variables
REAL, PARAMETER  :: ZSSTRENGTH = 0.3
REAL, PARAMETER  :: ZUSRMIN = 0.23
REAL, DIMENSION(SIZE(PUSR))  :: ZH
REAL, DIMENSION(SIZE(PUSR))  :: ZHH ! Gust level height (m)
REAL, DIMENSION(SIZE(PUSR))  :: ZRM
INTEGER :: JJ
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!---------------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:F09A',0,ZHOOK_HANDLE)
!
PSV(:) = 0.
PSA(:) = 0.
PFM(:) = 0. 
PVFM(:) = 0.
DO JJ=1,SIZE(PPHIOC)
   IF (PUSR(JJ) .GE. ZUSRMIN) THEN
!! and here, PHIOC = 0.1*ZU10 * ZUSR**2*ZRHOA/XRHOSW is expected
       ZH(JJ) =  PHS(JJ)/2.
       ZHH(JJ) = MAX(0.03*ZH(JJ),10.*PZOSEA(JJ))                                 
       PVFM(JJ) = (PVMOD(JJ) + PCP(JJ)/2.-PUSR(JJ)/0.4*LOG(PUREF(JJ)/ZHH(JJ))) * 0.07 * 1.15 - 0.3
       PVFM(JJ) = MAX(PVFM(JJ),0.04)
       ZRM(JJ) = 55.*PVFM(JJ)**0.7+20.
       PSV(JJ) = 1.E-5*(1.+(ZH(JJ)/3.)**0.1)*(ZRM(JJ)/50.)**2.5                    
       PSA(JJ) = 1.2 *(PPHIOC(JJ)/6E-4)**0.15*(ZRM(JJ)/73.)**(-1)  
       PFM(JJ) = XRHOSW*PSV(JJ)*PPHIOC(JJ)*ZSSTRENGTH
   ENDIF
END DO
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:F09A',1,ZHOOK_HANDLE)
END SUBROUTINE
!
!---------------------------------------------------------------------------------------
!#######################################################################################
SUBROUTINE F09C(PUSR,PU10,PUREF,PZOSEA,PHS,PCP,PRHOA,PSV,PSA,PVFM,PFM)
!#######################################################################################
!
!****  *F09C*
!
!       PURPOSE
!       -------
!       Computes key parameters of the sea-spray generation SUBROUTINE.
!       BASED ON BAO's FORTRAN CODE SHARED BY FAIRALL, V13
!       NB: PARAMS OF INTEGRATED SSGF STATS NOT GOING THROUGH ANY SIZE DISTRIBUTION
!         
!       Input Arguments:
!
!       Return:
!               PSV      --  normalized source volume (m/s)
!               PSA     --  normalized source area (1/s)
!               PVFM    --  estimate mean fall velocity (m/s)
!               PFM     --  droplet mass flux (kg/m2/s)
!
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       Fairall, C. W., Banner, M. L., Peirson, W. L., Asher, W., and Morison, R. P., 2009: 
!          Investigation of the physical scaling of sea spray spume droplet production, 
!          Journal of Geophysical Research, 114, C10 001, https://doi.org/10.1029/2008JC004918
!
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!
!       MODIFICATIONS
!       -------------
!---------------------------------------------------------------------------------------
!
USE MODD_OCEAN_CSTS, ONLY: XRHOSW
!
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
!
IMPLICIT NONE
!
!* 0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN) :: PUSR
REAL, DIMENSION(:), INTENT(IN) :: PU10
REAL, DIMENSION(:), INTENT(IN) :: PUREF
REAL, DIMENSION(:), INTENT(IN) :: PZOSEA
REAL, DIMENSION(:), INTENT(IN) :: PHS               ! Significant wave height (m)
REAL, DIMENSION(:), INTENT(IN) :: PCP               ! wave peak phase velocity (m/s)
REAL, DIMENSION(:), INTENT(IN) :: PRHOA
!
REAL, DIMENSION(:), INTENT(OUT) :: PFM    ! Spray mass flux (kg m-2 s-1)
REAL, DIMENSION(:), INTENT(OUT) :: PSV    ! normalized source volume m/s NB not used subsequently for now
REAL, DIMENSION(:), INTENT(OUT) :: PSA    ! normalized source area 1/s
REAL, DIMENSION(:), INTENT(OUT) :: PVFM   ! estimate mean fall velocity (m/s)
!
!* 0.2 declarations of local variables
REAL, PARAMETER :: ZSSTRENGTH = 0.3
REAL, DIMENSION(SIZE(PUSR))  :: ZH
REAL, DIMENSION(SIZE(PUSR))  :: ZHH ! Gust level height (m)
REAL, DIMENSION(SIZE(PUSR))  :: ZRM
REAL, DIMENSION(SIZE(PUSR))  :: ZPHIOC
REAL, PARAMETER :: ZU10MIN = 7.5
INTEGER :: JJ
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!---------------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:F09C',0,ZHOOK_HANDLE)
!
PSV(:) = 0.
PSA(:) = 0.
PFM(:) = 0.
PVFM(:) = 0.
DO JJ=1,SIZE(PUSR)
!! ad hoc ZPHIOC is used here
!!   ZPHIOC(JJ) = PRHOA(JJ)/XRHOSW*3.5*PUSR(JJ)**3.5  
   IF (PU10(JJ) .GE. ZU10MIN) THEN
       ZPHIOC(JJ) = PRHOA(JJ)*(-0.4+0.25*PU10(JJ))*PUSR(JJ)**2/XRHOSW
       ZH(JJ) =  PHS(JJ)/2.
       ZHH(JJ) = MAX(0.03*ZH(JJ),10.*PZOSEA(JJ))
       PVFM(JJ) = (1.2*PU10(JJ) - PCP(JJ)/2.-PUSR(JJ)/0.4*LOG(PUREF(JJ)/ZHH(JJ))) &
               *0.5*1.1/(1.+0.1*EXP(0.06*PU10(JJ)))
       PVFM(JJ) = MAX(PVFM(JJ),0.04)
       ZRM(JJ) = 183.*(PVFM(JJ)/1.32)**0.95
       PSV(JJ) = 5.E-5*(1.+(ZH(JJ)/3.)**0.1)*(ZRM(JJ)/183.)**2.5
       PSA(JJ) = 6.*(ZPHIOC(JJ)/6E-4)**0.3*(ZRM(JJ)/183.)**(-1)
       PFM(JJ) = ZSSTRENGTH*ZPHIOC(JJ)*XRHOSW*PSV(JJ)
   ENDIF
END DO
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:F09C',1,ZHOOK_HANDLE)
!
END SUBROUTINE
!
!---------------------------------------------------------------------------------------
!#######################################################################################
SUBROUTINE F94A(PU10, PRHOA, PVISA, PSV,PSA,PVFM, PFM)
!#######################################################################################
!
!****  *F94A*
!
!       PURPOSE
!       -------
!       Computes key parameters of the sea-spray generation SUBROUTINE.
!       Based on Fairall et al (1994)
!       Input Arguments:
!
!       Return:
!
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       Fairall, C. W., 1994
!          Investigation of the physical scaling of sea spray spume droplet production, 
!          Journal of Geophysical Research, 114, C10 001, https://doi.org/10.1029/2008JC004918
!
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!
!       MODIFICATIONS
!       -------------
!---------------------------------------------------------------------------------------
!
!USE MODD_SURF_PAR,       ONLY : XUNDEF
!
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
!
REAL, DIMENSION(:), INTENT(IN) :: PU10
REAL, DIMENSION(:), INTENT(IN) :: PVISA
REAL, DIMENSION(:), INTENT(IN) :: PRHOA
!
REAL, DIMENSION(:), INTENT(OUT) :: PFM    ! Spray mass flux (kg m-2 s-1)
REAL, DIMENSION(:), INTENT(OUT) :: PSV    ! normalized source volume m/s NB not used subsequently for now
REAL, DIMENSION(:), INTENT(OUT) :: PSA    ! normalized source area 1/s
REAL, DIMENSION(:), INTENT(OUT) :: PVFM   ! estimate mean fall velocity (m/s)
!
!* 0.2 declarations of local variables
!
INTEGER JJ, JI
INTEGER, PARAMETER :: ISIZE = 24
REAL  ::   ZRADIUS(ISIZE) = (/ 1., 2., 5., 10.,20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5 /) !in mum
REAL  ::   ZRADIUSM(ISIZE) = (/ 1., 2., 5., 10.,20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5 /)*1.E-6 !in m
REAL, DIMENSION(ISIZE)  :: ZR80
REAL, DIMENSION(ISIZE)  :: ZDR80DR0
REAL, PARAMETER  :: ZD1 = 1.02E4
REAL, PARAMETER  :: ZD2 = 6.95E6
REAL, PARAMETER  :: ZD3 = 1.75E17
REAL, PARAMETER  :: ZB0 = 4.405
REAL, PARAMETER  :: ZB1 = -2.646
REAL, PARAMETER  :: ZB2 = -3.156
REAL, PARAMETER  :: ZB3 = 8.902
REAL, PARAMETER  :: ZB4 = -4.482
REAL :: ZWINDFAC
REAL, PARAMETER :: ZU10MIN = 7.5
REAL, DIMENSION(ISIZE)  :: ZDFDR      ! Ejected droplet spectrum (SSGF) (m-2 s-1 m-1)
!---------------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:F94A',0,ZHOOK_HANDLE)
!
ZR80=0.518*ZRADIUS**0.976 !! Mueller & Veron 2014b eq 22
ZDR80DR0 = 0.506*ZRADIUS**(-0.024)
DO JJ = 1,SIZE(PU10)
    IF (PU10(JJ) .GE. ZU10MIN) THEN
        ZWINDFAC=(PU10(JJ)**3.4)/(11.**3.4) !! Mueller & Veron 2014b eq 20
        DO JI=1,ISIZE
           ZDFDR(JI) = ZWINDFAC*10**(ZB0 + ZB1*LOG10(ZR80(JI)) + ZB2*(LOG10(ZR80(JI)))**2 + &
                ZB3*(LOG10(ZR80(JI)))**3 + ZB4*(LOG10(ZR80(JI)))**4)
           IF (ZR80(JI) > 15.0) THEN
             ZDFDR(JI) = ZWINDFAC*ZD1/ZR80(JI)
             IF (ZR80(JI) > 37.0) THEN
                 ZDFDR(JI) = ZWINDFAC*ZD2*ZR80(JI)**(-2.8)
                 IF (ZR80(JJ) > 100.0) THEN
                     ZDFDR(JI) = ZWINDFAC*ZD3*ZR80(JI)**(-8.)
                 ENDIF
             ENDIF
           ENDIF
        ZDFDR(JI) = 1.E6*ZDR80DR0(JI) * ZDFDR(JI) !! M V eq 24
        END DO
!!! everything must be in m-2 s-1 m-1 (not mum-1) before calling SSGF_INT
        CALL SSGF_INT_STATS(ZDFDR, ZRADIUSM, PRHOA(JJ), PVISA(JJ), PSV(JJ), PSA(JJ), PVFM(JJ), PFM(JJ))
!!! consistency: droplets considered to be at equilibirum (R80), so integrate on r80 + seawater density
    ELSE 
        PVFM(JJ) = 0.
        PSV(JJ) = 0.
        PSA(JJ) = 0.
        PFM(JJ) = 0.
    ENDIF
END DO
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:F94A',1,ZHOOK_HANDLE)
!
END SUBROUTINE
!
!---------------------------------------------------------------------------------------
!#######################################################################################
SUBROUTINE B22B(PU10, PRHOA, PVISA, PSV,PSA,PVFM, PFM)
!#######################################################################################
!
!****  *B22B*
!
!       PURPOSE
!       -------
!       Computes key parameters of the sea-spray generation SUBROUTINE.
!       Based on Fairall et al (1994)
!       Modified as in Barr et al 2022
!       Input Arguments:
!
!       Return:
!
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       Fairall, C. W., 1994
!          Investigation of the physical scaling of sea spray spume droplet production, 
!          Journal of Geophysical Research, 114, C10 001, https://doi.org/10.1029/2008JC004918
!
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!
!       MODIFICATIONS
!       -------------
!
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
IMPLICIT NONE
!
!* 0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN) :: PU10
REAL, DIMENSION(:), INTENT(IN) :: PVISA
REAL, DIMENSION(:), INTENT(IN) :: PRHOA
!
REAL, DIMENSION(:), INTENT(OUT) :: PFM    ! Spray mass flux (kg m-2 s-1)
REAL, DIMENSION(:), INTENT(OUT) :: PSV    ! normalized source volume m/s NB not used subsequently for now
REAL, DIMENSION(:), INTENT(OUT) :: PSA    ! normalized source area 1/s
REAL, DIMENSION(:), INTENT(OUT) :: PVFM   ! estimate mean fall velocity (m/s)
!
!* 0.2 declarations of local variables
!
INTEGER JJ,JK
INTEGER, PARAMETER :: ISIZE = 24
!!REAL  ::   ZRADIUS(ISIZE) = (/ 1.,2., 5., 8., 10.,20.,25., 30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
!!        215.,250., 300.,350., 400.,450.,600. /) !in mum
!!REAL  ::   ZRADIUSM(ISIZE) = (/ 1.,2., 5.,8., 10.,20.,25.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
!!        215.,250., 300.,350., 400.,450.,600. /)*1.E-6 !in m
REAL  ::   ZRADIUS(ISIZE) = (/ 1., 2., 5., 10.,20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5 /) !in mum
REAL  ::   ZRADIUSM(ISIZE) = (/ 1., 2., 5., 10.,20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5 /)*1.E-6 !in m
REAL, DIMENSION(ISIZE)  ::   ZR80
REAL, DIMENSION(ISIZE)  ::   ZDR80DR0
REAL, PARAMETER :: ZD1 = 1.02E4
REAL, PARAMETER :: ZD2 = 6.95E6
REAL, PARAMETER :: ZD3 = 1.75E17
REAL, PARAMETER :: ZFS = 2.2
REAL, PARAMETER :: ZB0 = 4.405
REAL, PARAMETER :: ZB1 = -2.646
REAL, PARAMETER :: ZB2 = -3.156
REAL, PARAMETER :: ZB3 = 8.902
REAL, PARAMETER :: ZB4 = -4.482
REAL :: ZWINDFAC
REAL, PARAMETER :: ZU10MIN = 7.5
REAL :: ZWWI !! proportionality factor, as in Barr et al 2022
REAL, DIMENSION(ISIZE)  :: ZDFDR      ! Ejected droplet spectrum (SSGF) (m-2 s-1 m-1)
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!---------------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:B22B',0,ZHOOK_HANDLE)
!
ZR80=0.518*ZRADIUS**0.976 !! Mueller & Veron eq 22
ZDR80DR0 = 0.506*ZRADIUS**(-0.024)
DO JJ = 1,SIZE(PU10)
    IF (PU10(JJ) .GE. ZU10MIN) THEN
        ZWWI = (6.5E-4*(PU10(JJ)-2.)**1.5)/(3.8E-6*PU10(JJ)**3.4)
        ZWINDFAC=(PU10(JJ)**3.4)/(11.**3.4)
        DO JK=1,ISIZE
          ZDFDR(JK) = ZWINDFAC*10**(ZB0 + ZB1*LOG10(ZR80(JK)) + ZB2*(LOG10(ZR80(JK)))**2 + &
              ZB3*(LOG10(ZR80(JK)))**3 + ZB4*(LOG10(ZR80(JK)))**4)
          IF (ZR80(JK) > 15.0) THEN
             ZDFDR(JK) = ZWINDFAC*ZD1/ZR80(JK)
             IF (ZR80(JK) > 37.0) THEN
                 ZDFDR(JK) = ZWINDFAC*ZD2*ZR80(JK)**(-2.8)
                 IF (ZR80(JK) > 100.0) THEN
                     ZDFDR(JK) = ZWINDFAC*ZD3*ZR80(JK)**(-8.)
                 ENDIF
             ENDIF
          ENDIF
          ZDFDR(JK) = ZFS * ZWWI * ZDR80DR0(JK) * ZDFDR(JK) !! M V eq 24
        END DO
        ZDFDR = 1.E6*ZDFDR
!!! everything must be in m-2 s-1 m-1 (not mum-1) before calling SSGF_INT
        CALL SSGF_INT_STATS(ZDFDR, ZRADIUSM, PRHOA(JJ), PVISA(JJ), PSV(JJ), PSA(JJ), PVFM(JJ), PFM(JJ))
    ELSE 
        PVFM(JJ) = 0.
        PSV(JJ) = 0.
        PSA(JJ) = 0.
        PFM(JJ) = 0.
    ENDIF
END DO
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:B22B',1,ZHOOK_HANDLE)
END SUBROUTINE
!
!---------------------------------------------------------------------------------------
!#######################################################################################
SUBROUTINE F09B(PUSR, PU10, PZO, PHS, PCP, PMSS, PPHIOC, PRHOA, PVISA, PSV,PSA,PVFM, PFM)
!#######################################################################################
!
!****  *F09B*
!
!       PURPOSE
!       -------
!       Computes key parameters of the sea-spray generation SUBROUTINE.
!       Based on Fairall et al (2009), see Appendices, but directly taken from
!       the Matlab code of C. Fairall, version 7
!       Input Arguments:
!
!       Return:
!
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       Fairall, C. W., Banner, M. L., Peirson, W. L., Asher, W., and Morison, R. P., 2009: 
!          Investigation of the physical scaling of sea spray spume droplet production, 
!          Journal of Geophysical Research, 114, C10 001, https://doi.org/10.1029/2008JC004918
!
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!       MODIFICATIONS
!       -------------
!---------------------------------------------------------------------------------------
!
USE MODD_CSTS,       ONLY : XPI, XG, XKARMAN
!
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
!
IMPLICIT NONE
!
!* 0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN) :: PUSR
REAL, DIMENSION(:), INTENT(IN) :: PU10
REAL, DIMENSION(:), INTENT(IN) :: PHS               ! Significant wave height (m)
REAL, DIMENSION(:), INTENT(IN) :: PCP               ! wave peak phase velocity (m/s)
REAL, DIMENSION(:), INTENT(IN) :: PMSS              ! mean squared slope
REAL, DIMENSION(:), INTENT(IN) :: PPHIOC
REAL, DIMENSION(:), INTENT(IN) :: PVISA
REAL, DIMENSION(:), INTENT(IN) :: PRHOA
REAL, DIMENSION(:), INTENT(IN) :: PZO
!
REAL, DIMENSION(:), INTENT(OUT) :: PFM    ! Spray mass flux (kg m-2 s-1)
REAL, DIMENSION(:), INTENT(OUT) :: PSV    ! normalized source volume m/s NB not used subsequently for now
REAL, DIMENSION(:), INTENT(OUT) :: PSA    ! normalized source area 1/s
REAL, DIMENSION(:), INTENT(OUT) :: PVFM   ! estimate mean fall velocity (m/s)
!
!* 0.2 declarations of local variables
!
REAL  :: ZW          ! Active breaking whitecap fraction
REAL  :: ZH    ! 
REAL  :: ZHH    ! 
REAL  :: ZUD
REAL  :: ZDL
REAL  :: ZSIGUD
REAL  :: ZKP !!
INTEGER :: JJ
INTEGER, PARAMETER :: ISIZE = 25
REAL  ::   ZRADIUS(ISIZE) = (/ 10.,20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5,1250.,1500.,1750.,2000. /) !in mum
REAL  ::   ZRADIUSM(ISIZE) = (/ 10.,20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5,1250.,1500.,1750.,2000. /)*1.E-6 !in m
REAL, PARAMETER   ::   ZVISW  = 1.E-6                           ! Kinematic viscosity of water (m2 s-1)
REAL, PARAMETER   ::   ZSIGMA = 7.4E-5   ! Ratio of surface tension to water density [m3 s-2]
REAL, PARAMETER   ::   ZALPHA_K = 1.55    ! Kolmogorov constant [-], in F09
REAL, PARAMETER   ::   ZDLL = 1.0
REAL, PARAMETER   ::   ZFRAC = 0.12
REAL, PARAMETER   ::   ZU10MIN = 7.5
REAL  ::   ZND
REAL, DIMENSION(ISIZE)  :: ZVELF        ! Ejection probability (-)
REAL, DIMENSION(ISIZE)  :: ZDFDR      ! Ejected droplet spectrum (SSGF) (m-2 s-1 m-1)
REAL, DIMENSION(ISIZE)  :: ZVG      ! 
REAL, DIMENSION(ISIZE)  :: ZNLF
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!---------------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:F09B',0,ZHOOK_HANDLE)
!
!
DO JJ = 1,SIZE(PUSR)
   IF (PU10(JJ) .GE. ZU10MIN) THEN 
       ZW = MAX(2.E-7*PU10(JJ)**3.4, 1.E-3)
       ZH = PHS(JJ)/2.
       ZHH = MAX(0.03*ZH,10*PZO(JJ))
       ZKP = XG/PCP(JJ)**2
       ZUD = PUSR(JJ)/XKARMAN*LOG(ZHH/PZO(JJ))
       ZDL = 0.08 * (ZVISW**3/PPHIOC(JJ)*ZW)**0.25
!! here PHIOC is kept as it is, divided by XRHOSW in coare
       ZSIGUD = 0.75 * PU10(JJ)
       ZND = ZFRAC*3./(4.*XPI*ZSIGMA)
       CALL VG_PK97(ZRADIUSM, PRHOA(JJ), PVISA(JJ), ZVG) !  Settling velocity (m s-1)
!!   ZVELF = (ZUD + PCP(JJ)/2. - ZVG/PMSS(JJ)/1.3)/ZSIGUD    ! Ejection probability [-]
       ZVELF = (ZUD - PCP(JJ)/2. + 0.75*ZH*PCP(JJ) + 0.75*SQRT(XG/ZKP) - ZVG/PMSS(JJ)/1.3)/ZSIGUD    ! Ejection probability [-]
       ZNLF = (1.+ERF(ZVELF))/2.
!
       ZDFDR = ZND /ZRADIUSM**2*ZNLF*PPHIOC(JJ)/ZDLL*EXP(-9./4.*ZALPHA_K*(XPI*ZDL/ZRADIUSM)**(4./3.))
       CALL SSGF_INT_STATS(ZDFDR, ZRADIUSM, PRHOA(JJ), PVISA(JJ), PSV(JJ), PSA(JJ), PVFM(JJ), PFM(JJ))
       PFM(JJ) = MIN(PFM(JJ),1.E-1)
   ELSE
        PVFM(JJ) = 0.
        PFM(JJ) = 0.
        PSA(JJ) = 0.
        PSV(JJ) = 0.
   ENDIF
END DO
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:F09B',1,ZHOOK_HANDLE)
!
END SUBROUTINE
!
!---------------------------------------------------------------------------------------
!
!#######################################################################################
SUBROUTINE B22A(PUSR, PU10, PHS, PCP, PMSS, PPHIOC, PRHOA, PVISA, PSV,PSA,PVFM, PFM)
!#######################################################################################
!
!****  *B22A*
!
!       PURPOSE
!       -------
!       Computes key parameters of the sea-spray generation SUBROUTINE.
!       Based on Barr et al (2022).
!       Input Arguments:
!
!       Return:
!
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       Barr, B. W., Chen, S. S., & Fairall, C. W. (2022). Seastate-Dependent Sea Spray 
!             and Air-Sea Heat Fluxes in Tropical Cyclones: A New Parameterization for 
!             Fully Coupled Atmosphere-Wave-Ocean Models, 
!             Journal of the Atmospheric Sciences, DOI: 10.1175/JAS-D-22-0126.1
!             https://journals.ametsoc.org/view/journals/atsc/aop/JAS-D-22-0126.1/JAS-D-22-0126.1.xml
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!       MODIFICATIONS
!       -------------
!---------------------------------------------------------------------------------------
!
USE MODD_CSTS,       ONLY : XG, XPI, XKARMAN
USE MODD_OCEAN_CSTS, ONLY : XRHOSW
!
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
!
IMPLICIT NONE
!
!* 0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN) :: PUSR
REAL, DIMENSION(:), INTENT(IN) :: PU10
REAL, DIMENSION(:), INTENT(IN) :: PHS               ! Significant wave height (m)
REAL, DIMENSION(:), INTENT(IN) :: PCP               ! wave peak phase velocity (m/s)
REAL, DIMENSION(:), INTENT(IN) :: PMSS              ! mean squared slope
REAL, DIMENSION(:), INTENT(IN) :: PPHIOC
REAL, DIMENSION(:), INTENT(IN) :: PVISA
REAL, DIMENSION(:), INTENT(IN) :: PRHOA
!
REAL, DIMENSION(:), INTENT(OUT) :: PFM    ! Spray mass flux (kg m-2 s-1)
REAL, DIMENSION(:), INTENT(OUT) :: PSV    ! normalized source volume m/s NB not used subsequently for now
REAL, DIMENSION(:), INTENT(OUT) :: PSA    ! normalized source area 1/s
REAL, DIMENSION(:), INTENT(OUT) :: PVFM   ! estimate mean fall velocity (m/s)
!
!* 0.2 declarations of local variables
!
REAL, DIMENSION(SIZE(PUSR))  :: ZW          ! Active breaking whitecap fraction
REAL, DIMENSION(SIZE(PUSR))  :: ZPHIOC
REAL, DIMENSION(SIZE(PUSR))  :: ZEPS_KV     ! Surface vol kin dissipation under whitecaps (m2 s-3), 100 per S&M2015
REAL :: ZETA_K      ! Kolmogorov microscale (m)
REAL, DIMENSION(SIZE(PUSR))   :: ZU_H !  windspeed at gust height [m s-1]
REAL, DIMENSION(SIZE(PUSR))  :: ZU_CREST    ! Speed of breaking wave crests [m s-1], 0.8 per Banner et al. 2014
INTEGER, PARAMETER :: ISIZE=25
REAL  ::   ZRADIUS(ISIZE) = (/ 10.,20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5,1250.,1500.,1750.,2000. /) !in mum
REAL  ::   ZRADIUSM(ISIZE) = (/ 10.,20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5,1250.,1500.,1750.,2000. /)*1.E-6 !in m
REAL, PARAMETER   ::   ZVISW  = 1.E-6                           ! Kinematic viscosity of water (m2 s-1)
REAL, PARAMETER   ::   ZSIGMA = 7.4E-5   ! Ratio of surface tension to water density [m3 s-2]
REAL, PARAMETER   ::   ZALPHA_K = 1.5    ! Kolmogorov constant [-]REAL  ::   ZC1 = 1.35
REAL, PARAMETER   ::   ZSSTRENGHT=2.2
REAL, PARAMETER   ::   ZC1 = 1.35
REAL, PARAMETER   ::   ZC2 = 0.0518
REAL, PARAMETER   ::   ZC3 = 0.719
REAL, PARAMETER   ::   ZC4 = 2.17
REAL, PARAMETER   ::   ZC5 = 0.852
REAL, PARAMETER   ::   ZU10MIN = 7.5
REAL, DIMENSION(SIZE(ZRADIUS)) :: ZDMDRF     ! Spectrum of droplets formed (kg m-2 s-1 m-1)
REAL, DIMENSION(SIZE(ZRADIUS))  ::ZEP        ! Ejection probability (-)
REAL, DIMENSION(SIZE(ZRADIUS))  :: ZDMDR      ! Ejected droplet spectrum (SSGF) (kg m-2 s-1 m-1)
REAL, DIMENSION(SIZE(ZRADIUS))  :: ZDFDR      ! Ejected droplet spectrum (SSGF) (m-2 s-1 m-1)
REAL, DIMENSION(SIZE(ZRADIUS))  :: ZVG      ! 
INTEGER :: JJ
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!---------------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:B22A',0,ZHOOK_HANDLE)
!
ZPHIOC = XRHOSW*PPHIOC
ZU_CREST = 0.8*PCP    ! Speed of breaking wave crests [m s-1], 0.8 per Banner et al. 2014
ZW = MAX(1.E-3,MIN(0.18*PCP*PUSR**2/XG/PHS,1.0))  ! Active breaking whitecap fraction, per Deike et al. 2017
ZEPS_KV = 100.*ZPHIOC/PHS/XRHOSW/ZW    ! Surface vol kin dissipation under whitecaps [m2 s-3], 100 per S&M2015
ZU_H = PUSR/XKARMAN*LOG(200.)    ! Windspeed at gust height [m s-1]; assumes gust height of 200 x Z0
!
PVFM(:) = 0.
PFM(:) = 0.
PSA(:) = 0.
PSV(:) = 0.
DO JJ = 1,SIZE(PUSR)
    IF (PU10(JJ) .GE. ZU10MIN)  THEN
        CALL VG_PK97(ZRADIUSM, PRHOA(JJ), PVISA(JJ), ZVG) !  Settling velocity (m s-1)
        ZEP = (1. + ERF((ZU_H(JJ) - ZU_CREST(JJ) - ZVG/ZC3/PMSS(JJ))/ZC4/PU10(JJ) - ZC5))/2.    ! Ejection probability [-] -- U10 used instead of sigma_Uh
        ZETA_K = (ZVISW**3/ZEPS_KV(JJ))**0.25      ! Kolmogorov microscale [m]
        ZDMDRF = (ZSSTRENGHT*ZC1*XRHOSW*ZEPS_KV(JJ)*ZRADIUSM*ZW(JJ))/(3.*ZSIGMA) &
                 *EXP(-1.5*ZALPHA_K*ZC2*(XPI*ZETA_K/ZRADIUSM)**(4./3.))    ! Spectrum of droplets formed [kg m-2 s-1 m-1]    
        ZDMDR = ZDMDRF*ZEP    ! Ejected droplet spectrum (SSGF) (kg m-2 s-1 m-1)
        ZDFDR = ZDMDR/(XRHOSW*4./3.*XPI*ZRADIUSM**3)
        CALL SSGF_INT_STATS(ZDFDR, ZRADIUSM, PRHOA(JJ), PVISA(JJ), PSV(JJ), PSA(JJ), PVFM(JJ), PFM(JJ))
    ELSE
        PVFM(JJ) = 0.
        PFM(JJ) = 0.
        PSA(JJ) = 0.
        PSV(JJ) = 0.
    ENDIF
    PFM(JJ) = MIN(PFM(JJ),1.E-1)
END DO
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:B22A',1,ZHOOK_HANDLE)
!
END SUBROUTINE

!---------------------------------------------------------------------------------------
!
!#######################################################################################
SUBROUTINE B21A(PU10, PMSS, PVISA,PRHOA, PSV,PSA,PVFM, PFM)
!#######################################################################################
!
!****  *B21A*
!
!       PURPOSE
!       -------
!       Computes key parameters of the sea-spray generation SUBROUTINE.
!       Based on Bruch et al. (2021). and Bruch et al. 2023
!       Lognormal modes, scale parameter = MSS 
!       Input Arguments:
!
!       Return:
!
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!       MODIFICATIONS
!       -------------
!---------------------------------------------------------------------------------------
!
USE MODD_CSTS,       ONLY : XPI
!
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
IMPLICIT NONE
!
!* 0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN) :: PU10
REAL, DIMENSION(:), INTENT(IN) :: PMSS              ! mean squared slope
REAL, DIMENSION(:), INTENT(IN) :: PVISA             ! viscosity of air
REAL, DIMENSION(:), INTENT(IN) :: PRHOA
!
REAL, DIMENSION(:), INTENT(OUT) :: PFM    ! Spray mass flux (kg m-2 s-1)
REAL, DIMENSION(:), INTENT(OUT) :: PSV    ! normalized source volume m/s NB not used subsequently for now
REAL, DIMENSION(:), INTENT(OUT) :: PSA    ! normalized source area 1/s
REAL, DIMENSION(:), INTENT(OUT) :: PVFM   ! estimate mean fall velocity (m/s)
!
!* 0.2 declarations of local variables
INTEGER :: JI, JJ
INTEGER, PARAMETER :: ISIZE = 25
REAL  ::   ZRADIUS(ISIZE) = (/ 1., 2., 5., 10.,15., 20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5 /) * 0.2 !in mum
REAL  ::   ZRADIUSM(ISIZE) = (/ 1., 2., 5., 10.,15., 20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5 /)*0.2E-6 !in m
REAL, DIMENSION(ISIZE,3) :: ZSIGMA_I, ZF_I ,ZMU_I
REAL, DIMENSION(ISIZE)   :: ZR80
REAL, DIMENSION(ISIZE)   :: ZDR80DR0 
REAL, DIMENSION(ISIZE)   :: ZDFDR
REAL, PARAMETER   :: ZU10MIN = 7.5
!
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:B21A',0,ZHOOK_HANDLE)
!
ZR80=0.518*ZRADIUS**0.976 !! Mueller & Veron eq 22
ZDR80DR0 = 0.506*ZRADIUS**(-0.024)
ZSIGMA_I = SPREAD((/1.55,1.8,2.1 /),DIM=1 ,NCOPIES=ISIZE)
ZMU_I = SPREAD((/2.5,7.,25. /),DIM=1 ,NCOPIES=ISIZE)
ZF_I = SPREAD((/ 5.38E6*(PMSS)**2.45,1.94E6*(PMSS)**2.3,1.31E5*(PMSS)**2.39 /),DIM=1 ,NCOPIES=ISIZE)
DO JJ=1,SIZE(PMSS)
    IF (PU10(JJ) .GE. ZU10MIN) THEN
         ZDFDR=ZDR80DR0*SUM(ZF_I/((2.*XPI)**.5*LOG(ZSIGMA_I)) &
                     *EXP(-1./2.*(LOG(SPREAD(ZR80,DIM=2 ,NCOPIES=3)/ZMU_I)/LOG(ZSIGMA_I))**2),DIM=2)
         !   
         ZDFDR = 1.E6*ZDFDR
!!! everything must be in m-2 s-1 m-1 (not mum-1) before calling SSGF_INT
         CALL SSGF_INT_STATS(ZDFDR, ZRADIUSM, PRHOA(JJ), PVISA(JJ), PSV(JJ), PSA(JJ), PVFM(JJ), PFM(JJ))
    ELSE
        PVFM(JJ) = 0.
        PFM(JJ) = 0.
        PSA(JJ) = 0.
        PSV(JJ) = 0.
    ENDIF
END DO
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:B21A',1,ZHOOK_HANDLE)
!
END SUBROUTINE
!
!---------------------------------------------------------------------------------------
!
!#######################################################################################
SUBROUTINE B21B(PUSR, PMSS, PVISA,PRHOA, PSV,PSA,PVFM, PFM)
!#######################################################################################
!
!****  *B21B*
!
!       PURPOSE
!       -------
!       Computes key parameters of the sea-spray generation SUBROUTINE.
!       Based on Bruch et al. (2021). and Bruch et al 2023 
!       Lognormal modes, based on r80
!       Scaling parameter = Ps in Bruch
!       Input Arguments:
!
!       Return:
!
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!       MODIFICATIONS
!       -------------
!---------------------------------------------------------------------------------------
!
!
USE MODD_CSTS,       ONLY : XG, XPI
!
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
! 
IMPLICIT NONE
!
!* 0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN) :: PUSR
REAL, DIMENSION(:), INTENT(IN) :: PMSS             ! mean squared slope
REAL, DIMENSION(:), INTENT(IN) :: PVISA
REAL, DIMENSION(:), INTENT(IN) :: PRHOA
!
REAL, DIMENSION(:), INTENT(OUT) :: PFM    ! Spray mass flux (kg m-2 s-1)
REAL, DIMENSION(:), INTENT(OUT) :: PSV    ! normalized source volume m/s NB not used subsequently for now
REAL, DIMENSION(:), INTENT(OUT) :: PSA    ! normalized source area 1/s
REAL, DIMENSION(:), INTENT(OUT) :: PVFM   ! estimate mean fall velocity (m/s)
!
!* 0.2 declarations of local variables
INTEGER J,JI,JJ
INTEGER, PARAMETER :: ISIZE = 25
REAL  ::   ZRADIUS(ISIZE) = (/ 1., 2., 5., 10.,15., 20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5 /) * 0.2 !in mum
REAL  ::   ZRADIUSM(ISIZE) = (/ 1., 2., 5., 10.,15., 20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5 /)*0.2E-6 !in m
REAL, DIMENSION(ISIZE,3)  :: ZSIGMA_I, ZF_I ,ZMU_I
REAL, DIMENSION(ISIZE)   :: ZR80
REAL, DIMENSION(ISIZE)   :: ZDR80DR0
REAL, DIMENSION(ISIZE)  :: ZDFDR      ! Ejected droplet spectrum (SSGF) (m-2 s-1 m-1)
REAL  ::   ZPS
REAL, PARAMETER :: ZUSRMIN = 0.23
!
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:B21B',0,ZHOOK_HANDLE)
!
ZR80=0.518*ZRADIUS**0.976 !! Mueller & Veron eq 22
ZDR80DR0 = 0.506*ZRADIUS**(-0.024)
ZSIGMA_I = SPREAD((/1.55,1.8,2.1 /),DIM=1 ,NCOPIES=ISIZE)
ZMU_I = SPREAD((/2.5,7.,25. /),DIM=1 ,NCOPIES=ISIZE)
DO JJ=1,SIZE(PUSR)
    IF (PUSR(JJ) .GE. ZUSRMIN) THEN
         ZPS = PUSR(JJ)**3/XG/PVISA(JJ)*PMSS(JJ)
         IF (PUSR(JJ) > 1.2) ZPS = 1511.
!! hard limiter here, corresponds to 25 m/s, to avoid B21B to diverge by strong winds
         ZF_I = SPREAD((/ 4.76E1*ZPS**0.92,1.69*ZPS**1.41,4.5E-1*ZPS**1.11 /),DIM=1 ,NCOPIES=ISIZE)
         ZDFDR=ZDR80DR0*SUM(ZF_I/((2.*XPI)**.5*LOG(ZSIGMA_I)) &
                     *EXP(-1./2.*(LOG(SPREAD(ZR80,DIM=2 ,NCOPIES=3)/ZMU_I)/LOG(ZSIGMA_I))**2),DIM=2)
         !   
         ZDFDR = 1.E6*ZDFDR
!!! everything must be in m-2 s-1 m-1 (not mum-1) before calling SSGF_INT
         CALL SSGF_INT_STATS(ZDFDR, ZRADIUSM, PRHOA(JJ), PVISA(JJ), PSV(JJ), PSA(JJ), PVFM(JJ), PFM(JJ))
    ELSE
        PVFM(JJ) = 0.
        PFM(JJ) = 0.
        PSA(JJ) = 0.
        PSV(JJ) = 0.
    ENDIF
END DO
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:B21B',1,ZHOOK_HANDLE)
!
END SUBROUTINE
!
!---------------------------------------------------------------------------------------
!#######################################################################################
SUBROUTINE S20A(PU10,PVISA,PRHOA,PSV,PSA,PVFM, PFM)
!#######################################################################################
!
!****  *S20A*
!
!       PURPOSE
!       -------
!       Computes key parameters of the sea-spray generation SUBROUTINE.
!       Based on Shi et al. (2020). with whitecap fraction from Eq 17 
!       Input Arguments:
!
!       Return:
!
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       
!
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!       MODIFICATIONS
!       -------------
!---------------------------------------------------------------------------------------
!
!
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
IMPLICIT NONE
!
!* 0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN) :: PU10
REAL, DIMENSION(:), INTENT(IN) :: PVISA
REAL, DIMENSION(:), INTENT(IN) :: PRHOA
!
!
REAL, DIMENSION(:), INTENT(OUT) :: PFM    ! spray mass flux (kg m-2 s-1)
REAL, DIMENSION(:), INTENT(OUT) :: PSV    ! spray normalized source volume m/s NB not used subsequently for now
REAL, DIMENSION(:), INTENT(OUT) :: PSA    ! spray normalized source area 1/s
REAL, DIMENSION(:), INTENT(OUT) :: PVFM   ! spray mean fall velocity (m/s)
!
!* 0.2 declarations of local variables
!
INTEGER :: JJ, JI
INTEGER, PARAMETER :: ISIZE = 22
REAL  ::   ZRADIUS(ISIZE) = (/ 2., 5., 8., 10.,20.,25., 30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,250., 300.,350., 400.,450. /) !in mum
REAL  ::   ZRADIUSM(ISIZE) = (/ 2., 5.,8., 10.,20.,25.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,250., 300.,350., 400.,450. /)*1.E-6 !in m
REAL, DIMENSION(ISIZE)   :: ZDFDR
REAL  ::   ZW          ! Active breaking whitecap fraction
REAL  ::   ZA(4) = (/ 3.8E6, 9.5E3, 2.25E7, 2.08E19 /)
REAL  ::   ZB(4) = (/ -3., -1., -2.8, -8. /)
REAL  ::   ZRR(5) = (/ 2., 20., 75., 200., 500. /)
REAL, PARAMETER :: ZU10MIN = 7.5
!
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:S20A',0,ZHOOK_HANDLE)
!
DO JI = 1, SIZE(ZA)
    WHERE (ZRADIUS>=ZRR(JI) .AND. ZRADIUS<ZRR(JI+1))
        ZDFDR = ZA(JI)*ZRADIUS**ZB(JI)
    END WHERE
END DO
DO JJ = 1, SIZE(PU10)
    IF (PU10(JJ) .GE. 7.5) THEN
        ZW = MAX(1.E-3,MIN(3.84E-6*PU10(JJ)**3.41,1.0))    ! Eq 17 Shi 2020, cf Monahan 1980 eq 5
        ZDFDR = ZW*ZDFDR
        ZDFDR = 1.E6*ZDFDR
!!! everything must be in m-2 s-1 m-1 (not mum-1) before calling SSGF_INT
        CALL SSGF_INT_STATS(ZDFDR, ZRADIUSM, PRHOA(JJ), PVISA(JJ), PSV(JJ), PSA(JJ), PVFM(JJ), PFM(JJ))
    ELSE
        PVFM(JJ) = 0.
        PFM(JJ) = 0.
        PSA(JJ) = 0.
        PSV(JJ) = 0.
    ENDIF
END DO
!    
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:S20A',1,ZHOOK_HANDLE)
!
END SUBROUTINE
!
!---------------------------------------------------------------------------------------
!
!#######################################################################################
SUBROUTINE S20B(PU10,PVISA,PRHOA,PSV,PSA,PVFM, PFM)
!#######################################################################################
!
!****  *S20B*
!
!       PURPOSE
!       -------
!       Computes key parameters of the sea-spray generation SUBROUTINE.
!       Based on Shi et al. (2020). with whitecap fraction from Eq 18 
!       Input Arguments:
!
!       Return:
!
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       
!
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!       MODIFICATIONS
!       -------------
!---------------------------------------------------------------------------------------
!
!
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
IMPLICIT NONE
!
!* 0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN) :: PU10
REAL, DIMENSION(:), INTENT(IN) :: PVISA
REAL, DIMENSION(:), INTENT(IN) :: PRHOA
!
!
REAL, DIMENSION(:), INTENT(OUT) :: PFM    ! spray mass flux (kg m-2 s-1)
REAL, DIMENSION(:), INTENT(OUT) :: PSV    ! spray normalized source volume m/s NB not used subsequently for now
REAL, DIMENSION(:), INTENT(OUT) :: PSA    ! spray normalized source area 1/s
REAL, DIMENSION(:), INTENT(OUT) :: PVFM   ! spray mean fall velocity (m/s)
!
!* 0.2 declarations of local variables
!
INTEGER :: JJ, JI
INTEGER, PARAMETER :: ISIZE = 22
REAL  ::   ZRADIUS(ISIZE) = (/ 2., 5., 8., 10.,20.,25., 30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,250., 300.,350., 400.,450. /) !in mum
REAL  ::   ZRADIUSM(ISIZE) = (/ 2., 5.,8., 10.,20.,25.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,250., 300.,350., 400.,450. /)*1.E-6 !in m
REAL, DIMENSION(ISIZE)   :: ZDFDR
REAL  ::   ZW          ! Active breaking whitecap fraction
REAL  ::   ZA(4) = (/ 3.8E6, 9.5E3, 2.25E7, 2.08E19 /)
REAL  ::   ZB(4) = (/ -3., -1., -2.8, -8. /)
REAL  ::   ZRR(5) = (/ 2., 20., 75., 200., 500. /)
REAL, PARAMETER :: ZU10MIN = 7.5
!
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:S20B',0,ZHOOK_HANDLE)
!
DO JI = 1, SIZE(ZA)
    WHERE (ZRADIUS>=ZRR(JI) .AND. ZRADIUS<ZRR(JI+1))
        ZDFDR = ZA(JI)*ZRADIUS**ZB(JI)
    END WHERE
END DO
DO JJ = 1, SIZE(PU10)
    IF (PU10(JJ) .GE. ZU10MIN) THEN
        ZW = MAX(1.E-3,MIN(2.98E-5*PU10(JJ)**4.04,1.0))    ! Eq 18 Shi 2020, cf Zhao 2001 eq 30
        ZDFDR = ZW*ZDFDR
        ZDFDR = 1.E6*ZDFDR
!!! everything must be in m-2 s-1 m-1 (not mum-1) before calling SSGF_INT
        CALL SSGF_INT_STATS(ZDFDR, ZRADIUSM, PRHOA(JJ), PVISA(JJ), PSV(JJ), PSA(JJ), PVFM(JJ), PFM(JJ))
    ELSE
        PVFM(JJ) = 0.
        PFM(JJ) = 0.
        PSA(JJ) = 0.
        PSV(JJ) = 0.
    ENDIF
END DO
!    
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:S20B',1,ZHOOK_HANDLE)
!
END SUBROUTINE
!
!---------------------------------------------------------------------------------------
!#######################################################################################
SUBROUTINE T18A(PUSR,PVISA,PRHOA,PSV,PSA,PVFM, PFM)
!#######################################################################################
!
!****  *T18A*
!
!       PURPOSE
!       -------
!       Computes key parameters of the sea-spray generation SUBROUTINE.
!       Based on Troitskaya et al. (2018). 
!       Input Arguments:
!
!       Return:
!
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       Troitskaya, Y., A. Kandaurov, O. Ermakova, D. Kozlov, D. Sergeev, and S. Zilitinkevich, The Bag Breakup
!       Spume Droplet Generation Mechanism at High Winds. Part I: Spray Generation SUBROUTINE, Journal of Physical
!       Oceanography, 48 (9), 2167–2188, doi:10.1175/JPO-D-17-0104.1, 2018
!
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!       MODIFICATIONS
!       -------------
!---------------------------------------------------------------------------------------
!
USE MODD_CSTS,       ONLY : XG
USE MODD_OCEAN_CSTS, ONLY: XRHOSW
!
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
!
IMPLICIT NONE
!
!* 0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN) :: PUSR             ! friction velocity (m/s)
REAL, DIMENSION(:), INTENT(IN) :: PVISA            ! viscosity of air
REAL, DIMENSION(:), INTENT(IN) :: PRHOA            ! aie density (kg/m3)
!
!
REAL, DIMENSION(:), INTENT(OUT) :: PFM    ! spray mass flux (kg m-2 s-1)
REAL, DIMENSION(:), INTENT(OUT) :: PSV    ! spray normalized source volume m/s NB not used subsequently for now
REAL, DIMENSION(:), INTENT(OUT) :: PSA    ! spray normalized source area 1/s
REAL, DIMENSION(:), INTENT(OUT) :: PVFM   ! spray mean fall velocity (m/s)
!
!* 0.2 declarations of local variables
REAL    :: ZR2             ! radius at moment of rupture (m)
INTEGER :: JJ
INTEGER, PARAMETER :: ISIZE = 26
REAL  ::   ZRADIUS(ISIZE) = (/ 2.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,102.5,122.5,157.5,&
        215.,300.,400.,500.,600.,700.,800.,900.,1037.5,1250.,1500.,1800. /)*1.E-6 !in m
REAL, DIMENSION(ISIZE)   :: ZDFDR
REAL                     ::  ZTHETA1 !  ( )
REAL                     ::  ZTHETA2 ! (m)
REAL                     ::  ZN     ! number of bags ( )
REAL, PARAMETER          ::  ZSIGMA = 7.E-2 ! this value is for freshwater surface tension in N/m
REAL, PARAMETER          ::  ZL = 20. ! in [m]
REAL, PARAMETER          ::  ZNRIMD = 8.3   ! average number of rim droplets from one bag
REAL, PARAMETER          ::  ZUSRMIN = 0.23
REAL, PARAMETER  :: ZSSTRENGTH = 0.3
REAL, DIMENSION(SIZE(PUSR))  :: ZRM
!
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:T18A',0,ZHOOK_HANDLE)
!
DO JJ=1,SIZE(PUSR)
    IF (PUSR(JJ) .GE. ZUSRMIN) THEN
        ZR2 = 9.6/PUSR(JJ) * 1.E-3 ! EQ. 12, with ZR2 in mm, back to m here
        ZTHETA1 = 0.001*ZR2**(4./3.)/(ZL)**(1./3.) !! everything in m
        ZTHETA2 = 0.0021*ZR2 !! same
!!       ZN = 2.58E-4*ZRB**(3./2.)*EXP(-6.93E5/ZRB**(3./2.)) ! EQ. (7)
!!     MNB: WARNING EQ 7 not appropriate for field conditions; PTP is very different; revert to EQ 5
        ZN = 9.27E2*PUSR(JJ)**2.*EXP(-2.**2/PUSR(JJ)**2) ! EQ. (5)
        ZDFDR = ZN*(3.3E-9/(ZL)*(XRHOSW*XG*ZL**2/ZSIGMA)**(1.18)*(ZRADIUS/ZTHETA1)**7.3 &
                 *EXP(-5.2*SINH(3./7.*LOG(ZRADIUS/ZTHETA1)))                                &
                +1.5E-4*ZNRIMD/ZTHETA2*(ZRADIUS/ZTHETA2)**4.5                       &
                  *EXP(-3.94*SINH(1./2.*LOG(ZRADIUS/ZTHETA2)))) ! EQ. (24) everything in m 
!! already in m, not mum
        CALL SSGF_INT_STATS(ZDFDR, ZRADIUS, PRHOA(JJ), PVISA(JJ), PSV(JJ), PSA(JJ), PVFM(JJ), PFM(JJ))
    ELSE
        PVFM(JJ) = 0.
        PFM(JJ) = 0.
        PSA(JJ) = 0.
        PSV(JJ) = 0.
    ENDIF
END DO
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:T18A',1,ZHOOK_HANDLE)
!
END SUBROUTINE
!
!---------------------------------------------------------------------------------------
!
!#######################################################################################
SUBROUTINE VG_PK97(PRADIUS, PRHOA, PVISA, PVG)
!#######################################################################################
!
!****  *VG_PK97*
!
!       PURPOSE
!       -------
!       Fall velocity of spherical droplets.
!       Based on Pruppacher and Klett (1997) section 10.3.6.
!       Input Arguments:
!             PRADIUS - Droplet radius (m)
!             PRHOA   - 
!             PVISA   - Kinematic viscosity of air (m2 s-1)
!       Return:
!             PVG - settling velocity (m s-1)
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!       MODIFICATIONS
!       -------------
!-------------------------------------------------------------------------------
!
!USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODD_CSTS,       ONLY : XG
USE MODD_OCEAN_CSTS, ONLY: XRHOSW
!
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:),INTENT(IN) :: PRADIUS ! Droplet radius [m]
REAL, INTENT(IN) :: PRHOA
REAL, INTENT(IN) :: PVISA
!
REAL, DIMENSION(:),INTENT(OUT):: PVG
!
!* 0.2 declarations of local variables
!-------------------------------------------------------------------------------
!
REAL, PARAMETER :: ZSIGMA_AW = 7.4E-2    ! Surface tension of air-water interface [N m-1]
REAL, PARAMETER :: ZLAMDA_A0 = 6.6E-8    ! Mean free path at 1013.25 mb, 293.15 K [m]
REAL, DIMENSION(7) :: ZB1 = (/ -0.318657E1,0.992696,-0.153193E-2,-0.987059E-3,&
        -0.578878e-3,0.855176e-4,-0.327815e-5 /)    ! Polynomial coeffs for 10um < r < 535um
REAL, DIMENSION(6) :: ZB2 = (/ -0.500015E1,0.523778E1,-0.204914E1,0.475294,&
        -0.542819E-1,0.238449E-2 /)      ! Polynomial coeffs for r > 535um
REAL :: ZV_STOKES                        ! Stokes velocity [m s-1]
REAL :: ZF_SLIP                          ! Slip flow Cunningham correction factor [-]
REAL :: ZNRE1, ZNRE2                     ! Reynolds number [-]
REAL :: ZCDZNRE2                         ! Product of Cd and Re**2 [-]
REAL :: ZX1, ZX2                         ! 'X' in polynomial curve fit [-]
REAL :: ZY1, ZY2                         ! 'Y' in polynomial curve fit [-]
REAL :: ZNBO                             ! Bond number [-]
REAL :: ZNP                              ! Physical property number [-]
REAL :: ZNBONP16                         ! Product of ZNBO and ZNP**(1/6) [-]
INTEGER :: JJ
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:VG_PK97',0,ZHOOK_HANDLE)
!
DO JJ = 1, SIZE(PRADIUS)
    IF (PRADIUS(JJ) .LT. 10.E-6) THEN
        ZV_STOKES = 2.*PRADIUS(JJ)**2*XG*(XRHOSW - PRHOA)/9./(PRHOA*PVISA)    
        ZF_SLIP = 1. + 1.26*ZLAMDA_A0/PRADIUS(JJ)    
        PVG(JJ) = ZF_SLIP*ZV_STOKES    ! Settling velocity for r < 10um [m s-1]
    ELSE IF ((PRADIUS(JJ) .GE. 10.E-6) .AND. (PRADIUS(JJ) .LE. 535.E-6)) THEN
        ZCDZNRE2 = 32.*PRADIUS(JJ)**3*(XRHOSW - PRHOA)/PRHOA/PVISA**2*XG/3.    
        ZX1 = LOG(ZCDZNRE2)    
        ZY1 = ZB1(1) + ZB1(2)*ZX1 + ZB1(3)*ZX1**2 + ZB1(4)*ZX1**3 + &
             ZB1(5)*ZX1**4 + ZB1(6)*ZX1**5 + ZB1(7)*ZX1**6   
        ZNRE1 = EXP(ZY1)    
        PVG(JJ) = PVISA*ZNRE1/2./PRADIUS(JJ)    ! Settling velocity for 10 <= r <= 535 um [m s-1]
    ELSE IF (PRADIUS(JJ) .GT. 535.E-6) THEN
        ZNBO= XG*(XRHOSW - PRHOA)*PRADIUS(JJ)**2/ZSIGMA_AW    
        ZNP = ZSIGMA_AW**3/PRHOA**2/PVISA**4/XG/(XRHOSW - PRHOA)    
        ZNBONP16  = ZNBO*ZNP**(1./6.)    
        ZX2 = LOG(16./3.*ZNBONP16 )    
        ZY2 = ZB2(1) + ZB2(2)*ZX2 + ZB2(3)*ZX2**2 + ZB2(4)*ZX2**3 + &
             ZB2(5)*ZX2**4 + ZB2(6)*ZX2**5    
        ZNRE2 = ZNP**(1./6.)*EXP(ZY2)    
        PVG(JJ) = PVISA*ZNRE2/2./PRADIUS(JJ)    ! Settling velocity for r > 535 um [m s-1]
    END IF
END DO
!    
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:VG_PK97',1,ZHOOK_HANDLE)
!
END SUBROUTINE
!
!---------------------------------------------------------------------------------------
!
!#######################################################################################
SUBROUTINE SSGF_INT_STATS(PDFDR, PRADIUS, PRHOA, PVISA, PSV,PSA,PVFM, PFM)
!#######################################################################################
!
!****  **
!
!       PURPOSE
!       -------
!       Input Arguments:
!               PDFDR   --
!               PRADIUS -- droplet radius (m)
!               PRHOA   -- air density (kg/m3)
!               PVISA   -- air viscosity
!       Return:
!               PSV      --  normalized source volume (m/s)
!               PSA     --  normalized source area (1/s)
!               PVFM    --  estimate mean fall velocity (m/s)
!               PFM     --  droplet mass flux (kg/m2/s)
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       AUTHOR
!       ------
!       Sophia E. Brumer
!
!       MODIFICATIONS
!       -------------
!-------------------------------------------------------------------------------
!
!USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODD_CSTS,       ONLY : XG, XPI
USE MODD_OCEAN_CSTS, ONLY: XRHOSW
! 
!USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK
!USE PARKIND1, ONLY : JPRB
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:),INTENT(IN) :: PDFDR
REAL, DIMENSION(:),INTENT(IN) :: PRADIUS
REAL, INTENT(IN) :: PRHOA
REAL, INTENT(IN) :: PVISA     ! 
!
REAL, INTENT(OUT) :: PSV    ! normalized source volume (m/s) 
REAL, INTENT(OUT) :: PSA    ! normalized source area (1/s)
REAL, INTENT(OUT) :: PFM    ! droplet mass flux (kg/m2/s)
REAL, INTENT(OUT) :: PVFM   ! estimate mean fall velocity (m/s)
!
!* 0.2 declarations of local variables
!-------------------------------------------------------------------------------
!
REAL, DIMENSION(SIZE(PRADIUS)):: ZDELTAR
REAL, DIMENSION(SIZE(PRADIUS)):: ZVG
REAL :: ZMUL
INTEGER :: J
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:SSGF_INT_STATS',0,ZHOOK_HANDLE)
!
PSV = 0.
PSA = 0.
PFM = 0.
PVFM = 0.
!
J = SIZE(PRADIUS)
!! ZDELTAR should be in m, because PDFDR is in m-2 s-1 m-1
ZDELTAR(2:J-1) = .5*(PRADIUS(3:J)-PRADIUS(1:J-2))
ZDELTAR(1) = ZDELTAR(2)
ZDELTAR(J) = ZDELTAR(J-1)
PSV =   4./3.*XPI*DOT_PRODUCT((PRADIUS)**3*PDFDR,ZDELTAR)
PFM =  XRHOSW*PSV
!
PSA =  4.*XPI*DOT_PRODUCT((PRADIUS)**2*PDFDR,ZDELTAR) 
CALL VG_PK97(PRADIUS, PRHOA, PVISA, ZVG)
ZMUL=DOT_PRODUCT(PRADIUS**3*PDFDR,ZDELTAR)
IF (ZMUL .GT. 0) THEN
PVFM = DOT_PRODUCT(PRADIUS**3*PDFDR*ZVG,ZDELTAR)/ZMUL !! Everything is already in SI unts 
!!PVFM = DOT_PRODUCT((PRADIUS*1.E6)**3*PDFDR*ZVG,ZDELTAR)/DOT_PRODUCT((PRADIUS*1.E6)**3*PDFDR,ZDELTAR)
ENDIF
!    
!IF (LHOOK) CALL DR_HOOK('MODE_SSGF:SSGF_INT_STATS',1,ZHOOK_HANDLE)
!
END SUBROUTINE
!
END MODULE MODE_SSGF
