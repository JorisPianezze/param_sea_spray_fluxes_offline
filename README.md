This repository contains an offline interface to coupling_sltn.F90.

To use the code:
```bash
cd src/
make
./driver_param_sea_spray_fluxes_offline
```

This code reads a namelist called driver_param_sea_spray_fluxes_offline.nam :
```bash
&NAM_SURF_SLT CEMISPARAM_SLT = 'Ova14',
              LVARSIG_SLT    = .TRUE.,
              LRGFIX_SLT     = .FALSE.,
              JPMODE_SLT     = 5 /

&NAM_FORCING ZWIND  = 18.0,
             ZHS    = 2.0,
             ZSST   = 293.15,
             ZUSTAR = 1.6 /
```

| Namelist            | Variable name             |  Signification            |
| :------------------ | :------------------------ | :------------------------ |
| &NAM_SURF_SLT       | CEMISPARAM_SLT            | Type of parameterization  |
|                     | LVARSIG_SLT               | Sigma variable for of lognormal aerosol distributions |
|                     | LRGFIX_SLT                | Fix geometric radius of lognormal aerosol distributions  |
|                     | JPMODE_SLT                | Nb sea salt aerosol modes   |
| &NAM_FORCING        | ZWIND                     | Wind speed [m s-1]          |
|                     | ZHS                       | Significant wave height [m] |
|                     | ZSST                      | Sea surface temperature [K] |
|                     | ZUSTAR                    | Friction velocity [m s-1]   |
