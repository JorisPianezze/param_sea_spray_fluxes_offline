This repository contains an offline interface to coupling_sltn.F90.

To use the code:
```bash
cd src/
make
./driver_param_sea_spray_fluxes_offline
```

This code reads a file called driver_param_sea_spray_fluxes_offline.nam. This file contains two namelists described in the following table :


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
