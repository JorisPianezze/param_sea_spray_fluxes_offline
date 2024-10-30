!#########################################################
PROGRAM DRIVER_PARAM_SEA_SPRAY_FLUXES_OFFLINE
!#########################################################
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++ 
! ++   Interface offline to sea spray parametrizations :
! ++   - SURFEX/coupling_sltn.F90
! ++   - SURFEX/mode_ssgf.F90
! ++
! ++   original  : 07.04.2023 - J. Pianezze
! ++   revision  : 29.10.2024 - J. Pianezze
! ++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
USE MODD_SLT_n
USE MODD_SLT_SURF
USE MODD_CSTS,       ONLY : XG, XAVOGADRO, XPI, XKARMAN
USE MODD_OCEAN_CSTS, ONLY : XRHOSW 
USE MODD_SURFEX_n,   ONLY : SURFEX_t
USE MODN_SLT
USE MODE_SSGF
!
IMPLICIT NONE
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Sea spray & sea salt fluxes
!
REAL, DIMENSION(1  )              :: ZWIND       ! Wind speed [m s-1]
REAL, DIMENSION(1  )              :: ZHS         ! Significant wave height [m]
REAL, DIMENSION(1  )              :: ZUSTAR      ! Friction velocity [m s-1]
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Sea salt fluxes
!
INTEGER                           :: ISLT        ! Nb of sea salt moment variables
REAL, DIMENSION(1  )              :: ZSST        ! Sea surface temperature [K]
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSFSLT      ! Aerosol fluxes [kg m-2 s-1]
TYPE(SURFEX_t)                    :: YSURF_CUR   ! SURFEX derived data type   
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Sea spray fluxes
!
REAL, DIMENSION(1  )              :: ZCP          ! Wave peak phase velocity [m s-1]
REAL, DIMENSION(1  )              :: ZMSS         ! Mean squared slope
REAL, DIMENSION(1  )              :: ZPHIOC       ! ?
REAL, DIMENSION(1  )              :: ZRHOA        ! Air density [kg m-3]
REAL, DIMENSION(1  )              :: ZVISA        ! Air viscosity [kg m-3]
REAL, DIMENSION(1  )              :: ZSV          ! Normalized source volume [m s-1]  NB not used subsequently for now
REAL, DIMENSION(1  )              :: ZSA          ! Normalized source area [s-1]
REAL, DIMENSION(1  )              :: ZVFM         ! Estimate mean fall velocity [m s-1]
REAL, DIMENSION(1  )              :: ZFM          ! Spray mass flux [kg m-2 s-1]
!
NAMELIST /NAM_FORCING/ ZWIND, ZHS, ZSST, ZUSTAR 
!
! --------------------------------------------------------
!      1. Initialization
! --------------------------------------------------------
!
XG        = 9.80665 
XAVOGADRO = 6.0221367E+23
XPI       = 2.*ASIN(1.)
XKARMAN   = 0.4
XRHOSW    = 1024.
!
ZWIND (1) = 18.0
ZHS   (1) = 2.0
ZUSTAR(1) = 1.0
!
ZSST  (1) = 293.15
!
ZCP   (1) = 1.0
ZMSS  (1) = 0.3
ZPHIOC(1) = 0.2
ZRHOA (1) = 1.2
ZVISA (1) = 1E-5
!
CALL SLT_INIT(YSURF_CUR%SLT)
CALL INIT_SLT(YSURF_CUR%SLT,'MESONH')
!
! --------------------------------------------------------
!      2. Read namelists
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
CALL COUPLING_SLT_n (YSURF_CUR%SLT , &! [I  ] Sea salt (SLT) derived data type
                     SIZE(ZUSTAR,1), &! [I  ] Number of sea points 
                     ISLT,           &! [I  ] Number of sea salt variables 
                     ZWIND (:),      &! [I  ] Wind velocity [m s-1]
                     ZHS   (:),      &! [I  ] Significant wave height of wind-generated waves [m]
                     ZSST  (:),      &! [I  ] Sea surface temperature [K]
                     ZUSTAR(:),      &! [I  ] Friction velocity [m s-1]
                     ZSFSLT)          ! [  O] Production flux of sea salt [kg m-2 s-1] 
!
CLOSE(UNIT=11)
CLOSE(UNIT=12)
!
DEALLOCATE(ZSFSLT)
!
! --------------------------------------------------------
!      4. Compute sea spray aerosol fluxes
! --------------------------------------------------------
!
CALL B22A(ZUSTAR(:), &! [I  ] Friction velocity [m s-1] 
          ZWIND (:), &! [I  ] 10m wind speed [m s-1]
          ZHS   (:), &! [I  ] Significant wave height [m]
          ZCP   (:), &! [I  ] Wave peak phase velocity [m s-1]
          ZMSS  (:), &! [I  ] Mean squared slope
          ZPHIOC(:), &! [I  ] Wave to ocean energy flux
          ZRHOA (:), &! [I  ] Air density [kg m-3]
          ZVISA (:), &! [I  ] Air viscosity [kg m-3]
          ZSV   (:), &! [  O] Normalized source volume [m s-1] NB not used subsequently for now
          ZSA   (:), &! [  O] Normalized source area 1/s
          ZVFM  (:), &! [  O] Estimate mean fall velocity [m s-1]
          ZFM   (:))  ! [  O] Spray mass flux [kg m-2 s-1]
!
WRITE(*,*) ZSV, ZSA, ZVFM, ZFM
!
!#########################################################
END PROGRAM DRIVER_PARAM_SEA_SPRAY_FLUXES_OFFLINE
!#########################################################
