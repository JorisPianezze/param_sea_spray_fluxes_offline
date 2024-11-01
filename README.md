This repository contains an offline interface to coupling_sltn.F90 and mode_ssgf.F90.

# Python version

* Compilation :
  
```bash
cd src_python/SURFEX/
make
```

* Execution :
  
```bash
cd src_python/
python driver_param_sea_spray_fluxes_offline.py
```

or 

```bash
cd src_python/
jupyter-lab driver_param_sea_spray_fluxes_offline.ipynb
```
  
# Fortran version

* Compilation :
  
```bash
cd src_fortran/
make
```

* Execution :
  
```bash
cd src_fortran/
./driver_param_sea_spray_fluxes_offline
```

This executable reads a file called `driver_param_sea_spray_fluxes_offline.nam`. Namelist NAM_SURF_SLT is exactly the same as Meso-NH. This file contains two namelists described in the following table :

| Namelist      | Variable name  | Kind             | Signification                                                                                     |
|:------------- |:-------------- |:---------------- |:------------------------------------------------------------------------------------------------- |
| &NAM_SURF_SLT | CEMISPARAM_SLT | CHARACTER(LEN=6) | Type of parameterization ('Ova14', 'Vig01', 'OvB21a', 'OvB21b'). Default : 'Ova14'.               |
|               | LVARSIG_SLT    | LOGICAL          | Allow variation of sigma for of lognormal aerosol distributions. Default : .FALSE. (sigma fixed). |
|               | LRGFIX_SLT     | LOGICAL          | Fix geometric radius of lognormal aerosol distributions. Default : .FALSE. (radius variable).     |
|               | JPMODE_SLT     | INTEGER          | Nb of sea salt aerosol modes (3 to 8). Default : 5.                                               |
| &NAM_FORCING  | ZWIND          | REAL             | Wind speed [m s-1]. Default : 18.0.                                                               |
|               | ZHS            | REAL             | Significant wave height [m]. Default : 2.0.                                                       |
|               | ZSST           | REAL             | Sea surface temperature [K]. Default : 293.15                                                     |
|               | ZUSTAR         | REAL             | Friction velocity [m s-1]. Default : 1.6                                                          |

# Improvements

* ordre des variables dans b22a different de l'ordre d'appel dans la fonction F90 originale (pas un bug)

* ameliorer le jupyter notebook (commentaires, equations, etc...) pour le dossier example

* ajouter un fichier denvironnement conda avec les packages utilisés pour le developpement
