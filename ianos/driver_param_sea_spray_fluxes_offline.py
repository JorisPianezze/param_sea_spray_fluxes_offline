#!/usr/bin/python

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++
# ++   Interface offline to sea spray parametrizations :
# ++   - SURFEX/coupling_sltn.F90
# ++   - SURFEX/mode_ssgf.F90
# ++
# ++   original  : 07.04.2023 - J. Pianezze
# ++   revision  : 29.10.2024 - J. Pianezze
# ++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys, os
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from   matplotlib.colors import BoundaryNorm
from   mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
sys.path.append(os.path.abspath(os.curdir)+'/../src_python/SURFEX/')
import mode_ssgf
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# #########################################################
# ###           To be defined by user                   ###
# #########################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - -     File containing forcing fields                - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

cfg_file_WW3 = '/home/piaj/03_workdir/3M_study_IANOS/devel_SEASPRAY/ww3.202009.nc'

# #########################################################


# ---------------------------------------------------------
#      1. Initialization
# ---------------------------------------------------------

mode_ssgf.modd_csts.xg           = 9.80665
mode_ssgf.modd_csts.xpi          = 2*np.arcsin(1.0)
mode_ssgf.modd_csts.xkarman      = 0.4
mode_ssgf.modd_ocean_csts.xrhosw = 1024.0

# ---------------------------------------------------------
#      2. Read NetCDF file
# ---------------------------------------------------------

file_WW3 = netCDF4.Dataset(cfg_file_WW3)
lon_WW3  = file_WW3.variables['longitude']
lat_WW3  = file_WW3.variables['latitude']
nlon_WW3 = np.size(file_WW3.dimensions['longitude'])
nlat_WW3 = np.size(file_WW3.dimensions['latitude'])

wind_WW3  = file_WW3.variables['hs'] [80,:,:].reshape(nlon_WW3*nlat_WW3)+10.0 # BUG
hs_WW3    = file_WW3.variables['hs'] [80,:,:].reshape(nlon_WW3*nlat_WW3)
ustar_WW3 = file_WW3.variables['hs'] [80,:,:].reshape(nlon_WW3*nlat_WW3) # BUG
cp_WW3    = file_WW3.variables['tp'] [80,:,:].reshape(nlon_WW3*nlat_WW3)*(9.81/(2.0*np.pi))
mss_WW3   = file_WW3.variables['hs'] [80,:,:].reshape(nlon_WW3*nlat_WW3) # BUG
phioc_WW3 = file_WW3.variables['foc'][80,:,:].reshape(nlon_WW3*nlat_WW3)
rhoa_WW3  = file_WW3.variables['hs'] [80,:,:].reshape(nlon_WW3*nlat_WW3) # BUG
visa_WW3  = file_WW3.variables['hs'] [80,:,:].reshape(nlon_WW3*nlat_WW3)*1E-5 # BUG

# ---------------------------------------------------------
#      3. Compute sea spray aerosol fluxes
# ---------------------------------------------------------

sv_WW3,sa_WW3,vfm_WW3,fm_WW3 = mode_ssgf.mode_ssgf.b22a(ustar_WW3,wind_WW3,hs_WW3,cp_WW3,mss_WW3,phioc_WW3,rhoa_WW3,visa_WW3)

# ---------------------------------------------------------
#      4. Reshape variables
# ---------------------------------------------------------

wind_WW3  = wind_WW3.reshape (nlat_WW3,nlon_WW3)
hs_WW3    = hs_WW3.reshape   (nlat_WW3,nlon_WW3)
ustar_WW3 = ustar_WW3.reshape(nlat_WW3,nlon_WW3)
cp_WW3    = cp_WW3.reshape   (nlat_WW3,nlon_WW3)
mss_WW3   = mss_WW3.reshape  (nlat_WW3,nlon_WW3)
phioc_WW3 = phioc_WW3.reshape(nlat_WW3,nlon_WW3)
rhoa_WW3  = rhoa_WW3.reshape (nlat_WW3,nlon_WW3)
visa_WW3  = visa_WW3.reshape (nlat_WW3,nlon_WW3)

sv_WW3    = sv_WW3.reshape   (nlat_WW3,nlon_WW3)
sa_WW3    = sa_WW3.reshape   (nlat_WW3,nlon_WW3)
vfm_WW3   = vfm_WW3.reshape  (nlat_WW3,nlon_WW3)
fm_WW3    = fm_WW3.reshape   (nlat_WW3,nlon_WW3)

# ---------------------------------------------------------
#      5. Figure
# ---------------------------------------------------------

fig = plt.figure(figsize = (6,5))
gs  = gridspec.GridSpec(1,2)

# ---------------------------------------------------------
#    Add basemap
# ---------------------------------------------------------

bmap=Basemap(llcrnrlat=np.min(lat_WW3),llcrnrlon=np.min(lon_WW3),\
             urcrnrlat=np.max(lat_WW3),urcrnrlon=np.max(lon_WW3),\
             lat_1=np.min(lat_WW3),lat_2=np.max(lat_WW3),\
             lat_0=(np.min(lat_WW3)+np.max(lat_WW3))/2.0,\
             lon_0=(np.min(lon_WW3)+np.max(lon_WW3))/2.0,\
             projection='lcc',suppress_ticks=True, resolution = 'i')
    
mer=np.arange(-150, 200, 5.0)
par=np.arange(-120, 200, 5.0)
    
lon_WW3, lat_WW3 = bmap(lon_WW3, lat_WW3)

# ---------------------------------------------------------
#           ___________                  
#          |     |     |                  
#          |  X  |     |                  
#          |_____|_____|                  
#                                         
# ---------------------------------------------------------

ax00  = fig.add_subplot(gs[0,0])

# ---------------------------------------------------------
#    Add basemap info
# ---------------------------------------------------------

bmap.drawmeridians(mer, labels=[0, 0, 0, 1], fontsize=8, color='0.6', linewidth=0.5, ax=ax00)
bmap.drawparallels(par, labels=[1, 0, 0, 0], fontsize=8, color='0.6', linewidth=0.5, ax=ax00)
bmap.drawcoastlines(linewidth=0.5, linestyle='solid', color='k',ax=ax00)

# ---------------------------------------------------------
#    Definition de la colormap
# ---------------------------------------------------------

cmap_hs=plt.get_cmap('hot_r')
    
# ---------------------------------------------------------
#    Definition des niveaux
# ---------------------------------------------------------

levels_hs = np.arange(0.0, 5.1, 0.1)
norm_hs   = BoundaryNorm(levels_hs, ncolors=cmap_hs.N, clip=True)
    
# ---------------------------------------------------------
#    Pcolormesh
# ---------------------------------------------------------

bmap_pmsh = bmap.pcolormesh(lon_WW3[:,:],lat_WW3[:,:],hs_WW3[:,:],cmap=cmap_hs, norm=norm_hs)
cbar      = bmap.colorbar(bmap_pmsh) 

# ---------------------------------------------------------
#    Add texts
# ---------------------------------------------------------

plt.text(0.0,1.0,r'(a) Hs [m]',\
         horizontalalignment='left',verticalalignment='bottom',\
         fontsize=8, transform = ax00.transAxes, color='k')

# ---------------------------------------------------------
#           ___________                  
#          |     |     |                  
#          |     |  X  |                  
#          |_____|_____|                  
#                                         
# ---------------------------------------------------------


ax01  = fig.add_subplot(gs[0,1])

# ---------------------------------------------------------
#    Add basemap info
# ---------------------------------------------------------

bmap.drawmeridians(mer, labels=[0, 0, 0, 1], fontsize=8, color='0.6', linewidth=0.5, ax=ax01)
bmap.drawparallels(par, labels=[0, 0, 0, 0], fontsize=8, color='0.6', linewidth=0.5, ax=ax01)
bmap.drawcoastlines(linewidth=0.5, linestyle='solid', color='k', ax=ax01)

# ---------------------------------------------------------
#    Definition de la colormap
# ---------------------------------------------------------

cmap_flux=plt.get_cmap('hot_r')
    
# ---------------------------------------------------------
#    Definition des niveaux
# ---------------------------------------------------------

levels_flux = np.arange(0.0, 0.1, 0.01)
norm_flux   = BoundaryNorm(levels_flux, ncolors=cmap_flux.N, clip=True)
    
# ---------------------------------------------------------
#    Pcolormesh
# ---------------------------------------------------------

bmap_pmsh = bmap.pcolormesh(lon_WW3[:,:], lat_WW3[:,:], fm_WW3[:,:], cmap=cmap_flux, norm=norm_flux)
cbar      = bmap.colorbar(bmap_pmsh) 

# ---------------------------------------------------------
#    Add texts
# ---------------------------------------------------------

plt.text(0.0,1.0,r'(b) Mass flux [kg m-2 s-1]',\
         horizontalalignment='left',verticalalignment='bottom',\
         fontsize=8, transform = ax01.transAxes, color='k')

# ---------------------------------------------------------
#    Save figure
# ---------------------------------------------------------

plt.savefig('hs_mass_flux.png', bbox_inches='tight', dpi=400)
plt.close()
