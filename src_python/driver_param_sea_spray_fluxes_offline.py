#!/usr/bin/python

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++
# ++   Interface offline to sea spray parametrizations :
# ++   - SURFEX/coupling_sltn.F90
# ++   - SURFEX/mode_ssgf.F90
# ++
# ++   original  : 07.04.2023 - J. Pianezze
# ++   revision  : 29.10.2024 - J. Pianezze
# ++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys, os
import numpy as np
sys.path.append(os.path.abspath(os.curdir)+'/SURFEX/')
import mode_ssgf
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --------------------------------------------------------
#      1. Initialization
# --------------------------------------------------------

mode_ssgf.modd_csts.xg           = 9.80665
mode_ssgf.modd_csts.xpi          = 2*np.arcsin(1.0)
mode_ssgf.modd_csts.xkarman      = 0.4
mode_ssgf.modd_ocean_csts.xrhosw = 1024.0

wind  = 18.0
hs    = 2.0
ustar = 1.6
cp    = 1.0
mss   = 0.3
phioc = 0.2
rhoa  = 1.2
visa  = 1E-5

# --------------------------------------------------------
#      2. Compute sea salt aerosol fluxes
# --------------------------------------------------------

# ...

# --------------------------------------------------------
#      3. Compute sea spray aerosol fluxes
# --------------------------------------------------------

sv,sa,vfm,fm = mode_ssgf.mode_ssgf.b22a(ustar,wind,hs,cp,mss,phioc,rhoa,visa)

print(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
print(' ~~~                Sea spray fluxes                        ~~~ ')
print('validation                  : ', sv[0],sa[0],vfm[0],fm[0])
print('validation from src_fortran : ', 1.28740892E-06,67.3112488,6.57372856,0.100000001)
print(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
