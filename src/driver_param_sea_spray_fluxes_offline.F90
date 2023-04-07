!#########################################################
PROGRAM DRIVER_PARAM_BULK_OWA_OFFLINE
!#########################################################
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++ 
! ++   Programme permettant de tester les parametrisations
! ++   des flux a l'interface air/mer en offline
! ++
! ++   original  : 02.04.2021 - J. Pianezze
! ++   revision  : XX.XX.XXXX - XXXXXXXXXXX
! ++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
USE MODD_SLT_n
USE MODD_SLT_SURF
USE MODD_CSTS, ONLY : XAVOGADRO, XPI
USE MODD_SURFEX_n, ONLY : SURFEX_t
!
IMPLICIT NONE
!
! --------------------------------------------------------
! --
! --   0. Declaration des variables
! --  
! --------------------------------------------------------
!
! --   ---------------------------------
! --   0.1  Entrees des parametrisations
! --   ---------------------------------
!
REAL, DIMENSION(1) :: XTA       ! air temperature at atm. level (K)
REAL, DIMENSION(1) :: XQA       ! air humidity at atm. level (kg/kg)
REAL, DIMENSION(1) :: XPA       ! air pressure at atm level (Pa)
REAL, DIMENSION(1) :: XEXNA     ! Exner function at atm. level
REAL, DIMENSION(1) :: XRHOA     ! air density at atm. level
REAL, DIMENSION(1) :: XVMOD     ! module of wind at atm. wind level (m/s)
REAL, DIMENSION(1) :: XZREF     ! atm. level for temp. and humidity (m)
REAL, DIMENSION(1) :: XUREF     ! atm. level for wind (m)
REAL, DIMENSION(1) :: XSST      ! Sea Surface Temperature (K)
REAL, DIMENSION(1) :: XSSS      ! Sea Surface Salinity (g/kg)
REAL, DIMENSION(1) :: XEXNS     ! Exner function at sea surface
REAL, DIMENSION(1) :: XPS       ! air pressure at sea surface (Pa)
REAL, DIMENSION(1) :: XRAIN     ! precipitation rate (kg/s/m2)
REAL, DIMENSION(1) :: XPERTFLUX ! stochastic flux perturbation pattern
REAL               :: XICHCE    !
LOGICAL            :: LPRECIP   !
LOGICAL            :: LPWEBB    !
LOGICAL            :: LPWG      !
LOGICAL            :: LPERTFLUX !
INTEGER            :: NZ0       !
!
! --   ---------------------------------
! --   0.2  Sorties des parametrisations
! --   ---------------------------------
!
REAL, DIMENSION(1) :: XSFTH   ! heat flux (W/m2)
REAL, DIMENSION(1) :: XSFTQ   ! water flux (kg/m2/s)
REAL, DIMENSION(1) :: XUSTAR  ! friction velocity (m/s)
REAL, DIMENSION(1) :: XRI,XRESA,XQSAT,XCH,XCE,XCD,XCDN
REAL, DIMENSION(1) :: XZ0SEA  ! dynamical roughness length over the ocean (m)
REAL, DIMENSION(1) :: XZ0HSEA ! thermical roughness length over the ocean (m)
REAL, DIMENSION(1) :: XHS
REAL, DIMENSION(1) :: XTP
!
REAL, DIMENSION(1,23) :: XSFSLT      !Out: kg/m2/s (index #2)
!
! --   ---------------------------------------------
! --   0.3  Constantes, indice de boucle et namelist
! --   ---------------------------------------------
!
INTEGER            :: JWIND, JITER, NITERMAX
!
TYPE(SURFEX_t) :: YSURF_CUR
!
NAMELIST /NAM_BULK_FLUX/ LPRECIP, LPWEBB, LPWG, LPERTFLUX, &
                         XPS, XTA, XSST, XQA, XSSS
!
! --------------------------------------------------------
! --
! --   1. Initialisation des variables
! --  
! --------------------------------------------------------
!
XAVOGADRO   = 6.0221367E+23
XPI         = 2.*ASIN(1.)
!
LPRECIP   = .FALSE.
LPWEBB    = .FALSE.
LPWG      = .FALSE.
LPERTFLUX = .FALSE.
XICHCE    = 0.0
NZ0       = 0
XHS       = 0.0
XTP       = 0.0

XPS(1)    = 100000.

XRAIN(1)  = 0.02

XZREF(1)  = 10.
XUREF(1)  = 10.

XTA(1)    = 26. + 273.16
XQA(1)    = 20.E-3
XSST(1)   = 28. + 273.16
XSSS(1)   = 35.0

! --------------------------------------------------------
! --
! --   2. Lecture de la namelist
! --  
! --------------------------------------------------------
!
OPEN(UNIT=10,FILE='driver_param_sea_spray_fluxes_offline.nam')
READ(UNIT=10,NML=NAM_BULK_FLUX)
CLOSE(10)
                    
! --------------------------------------------------------
! --
! --   4. Calcul des flux pour les differentes param.
! --  
! --------------------------------------------------------
!
!OPEN(UNIT=12,FILE="SFTH")
!OPEN(UNIT=13,FILE="SFTQ")
!OPEN(UNIT=14,FILE="TAU")
OPEN(UNIT=12,FILE="XSFTSLT")
!
XVMOD(1) = 0.0
!
CEMISPARAM_SLT="Ova14"
JPMODE_SLT=8
!
CALL SLT_INIT(YSURF_CUR%SLT)
CALL INIT_SLT(YSURF_CUR%SLT,'MESONH')
!
!DO JWIND = 1, 61
DO JWIND = 1, 1
  !
  print*,' ~~ dans driver_bulk_flux.f90 2.0 : XVMOD=',XVMOD(1)
  !
  CALL COUPLING_SLT_n (YSURF_CUR%SLT, &
      1,                       &!I [nbr] number of sea points 
      23,                        &!I [nbr] number of sea salt variables 
      10.0,                    &!I Wind velocity
      2.0,                 &! Significant height of wind-generated waves (in ECMWF analyses)                     ! local pour l'instant, PWHEIGHT plus tard
      293.15,                     &! Sea water temperature (C) 
      0.2,                   &! Friction velocity (ecmwf?) Calcule dans coupling_seafluxn.F90
      XSFSLT)                    !O [kg/m2/sec] production flux of sea salt
  !
  XVMOD(1) = XVMOD(1) + 1.
  !
ENDDO
!
CLOSE(UNIT=12)
!
!#########################################################
END PROGRAM DRIVER_PARAM_BULK_OWA_OFFLINE
!#########################################################
