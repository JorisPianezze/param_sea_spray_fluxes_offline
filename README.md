This repository contains an offline interface to coupling_sltn.F90.

To use the code:
```bash
cd src/
make
./driver_param_sea_spray_fluxes_offline
```

This code reads a file called driver_param_sea_spray_fluxes_offline.nam. This file contains two namelists described in the following table :


| Namelist            | Name             |  Kind | Signification            |
| :------------------ | :------------------------ | :------------ |  :------------------------ |
| &NAM_SURF_SLT       | CEMISPARAM_SLT            | Character(\*8) | Type of parameterization  |
|                     | LVARSIG_SLT               | Logical | Sigma variable for of lognormal aerosol distributions |
|                     | LRGFIX_SLT                | Logical | Fix geometric radius of lognormal aerosol distributions  |
|                     | JPMODE_SLT                | Integer | Nb sea salt aerosol modes   |
| &NAM_FORCING        | ZWIND                     | Real | Wind speed [m s-1]          |
|                     | ZHS                       | Real | Significant wave height [m] |
|                     | ZSST                      | Real | Sea surface temperature [K] |
|                     | ZUSTAR                    | Real | Friction velocity [m s-1]   |
