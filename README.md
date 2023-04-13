This repository contains an offline interface to coupling_sltn.F90.

To use the code:
```bash
cd src/
make
./driver_param_sea_spray_fluxes_offline
```

This code reads a file called driver_param_sea_spray_fluxes_offline.nam. This file contains two namelists described in the following table :


| Namelist            | Variable name             |  Kind            | Signification            |
| :------------------ | :------------------------ | :--------------- |  :------------------------ |
| &NAM_SURF_SLT       | CEMISPARAM_SLT            | CHARACTER(LEN=6) | Type of parameterization  |
|                     | LVARSIG_SLT               | LOGICAL          | Sigma variable for of lognormal aerosol distributions |
|                     | LRGFIX_SLT                | LOGICAL          | Fix geometric radius of lognormal aerosol distributions  |
|                     | JPMODE_SLT                | INTEGER          | Nb of sea salt aerosol modes   |
| &NAM_FORCING        | ZWIND                     | REAL             | Wind speed [m s-1]          |
|                     | ZHS                       | REAL             | Significant wave height [m] |
|                     | ZSST                      | REAL             | Sea surface temperature [K] |
|                     | ZUSTAR                    | REAL             | Friction velocity [m s-1]   |
