!#########################################################
PROGRAM DRIVER_PARAM_SEA_SPRAY_FLUXES_OFFLINE
!#########################################################
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++ 
! ++   Interface offline to sea salt aerosol parametrization
! ++   SURFEX/coupling_sltn.F90
! ++
! ++   original  : 07.04.2023 - J. Pianezze
! ++   revision  : XX.XX.XXXX - XXXXXXXXXXX
! ++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
USE MODD_SLT_n
USE MODD_SLT_SURF
USE MODD_CSTS, ONLY : XAVOGADRO, XPI
USE MODD_SURFEX_n, ONLY : SURFEX_t
USE MODN_SLT
!
IMPLICIT NONE
!
INTEGER                           :: ISLT        ! Nb of sea salt moment variables
REAL, DIMENSION(1  )              :: ZWIND       ! Wind speed [m s-1]
REAL, DIMENSION(1  )              :: ZHS         ! Significant wave height [m]
REAL, DIMENSION(1  )              :: ZSST        ! Sea surface temperature [K]
REAL, DIMENSION(1  )              :: ZUSTAR      ! Friction velocity [m s-1]
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSFSLT      ! Aerosol fluxes [kg m-2 s-1]
TYPE(SURFEX_t)                    :: YSURF_CUR   ! SURFEX derived data type                   
!
NAMELIST /NAM_FORCING/ ZWIND, ZHS, ZSST, ZUSTAR 
!
! --------------------------------------------------------
!      1. Initialization
! --------------------------------------------------------
!
XAVOGADRO = 6.0221367E+23
XPI       = 2.*ASIN(1.)
ZWIND (1) = 18.0
ZHS   (1) = 2.0
ZSST  (1) = 293.15
ZUSTAR(1) = 1.0
!
CALL SLT_INIT(YSURF_CUR%SLT)
CALL INIT_SLT(YSURF_CUR%SLT,'MESONH')
!
! --------------------------------------------------------
!      2. Read namelist NAM_SURF_SLT
! --------------------------------------------------------
!
OPEN(UNIT=10, FILE='driver_param_sea_spray_fluxes_offline.nam')
READ(UNIT=10, NML=NAM_SURF_SLT)
READ(UNIT=10, NML=NAM_FORCING)
CLOSE(10)
!
IF (LVARSIG_SLT) THEN     ! 3 moments for each mdoe
  ISLT = JPMODE_SLT*3
ELSEIF (LRGFIX_SLT) THEN  ! 1 moment  for each mode
  ISLT = JPMODE_SLT
ELSE                      ! 2 moments for each mode
  ISLT = JPMODE_SLT*2
ENDIF
!
! --------------------------------------------------------
!      3. Compute sea salt aerosol fluxes
! --------------------------------------------------------
!
ALLOCATE(ZSFSLT(1,ISLT))
!
OPEN(UNIT=11, FILE="ZSFSLT_MDE")
OPEN(UNIT=12, FILE="ZSFSLT")
!
!
CALL COUPLING_SLT_n (YSURF_CUR%SLT , &
                     SIZE(ZUSTAR,1), &! [I  ] Number of sea points 
                     ISLT,           &! [I  ] Number of sea salt variables 
                     ZWIND(1),       &! [I  ] Wind velocity [m s-1]
                     ZHS(1),         &! [I  ] Significant wave height of wind-generated waves [m] ?
                     ZSST(1),        &! [I  ] Sea surface temperature [K]
                     ZUSTAR(1),      &! [I  ] Friction velocity (ecmwf?) Calcule dans coupling_seafluxn.F90 ?
                     ZSFSLT)          ! [  O] Production flux of sea salt [kg/m2/sec] 
!
CLOSE(UNIT=11)
CLOSE(UNIT=12)
!
DEALLOCATE(ZSFSLT)
!
!#########################################################
END PROGRAM DRIVER_PARAM_SEA_SPRAY_FLUXES_OFFLINE
!#########################################################
