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

|               &NAM_SURF_SLT                   ||
|                  :----:                       ||
| Variable            |  Signification           |
| :---                |                     ---: |
| CEMISPARAM_SLT      | Type of parameterization |
