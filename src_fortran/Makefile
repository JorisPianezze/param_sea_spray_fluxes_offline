##########################################################
#
# To compile the code :
# make
#
# To clean the compile directory :
# make clean
#
##########################################################

FFLAGS=-O0 
#-fimplicit-none -Wall
CPPFLAGS=-I${PWD}

default:
	gfortran -c -g ${FFLAGS} -o SURFEX/modd_csts.o                      SURFEX/modd_csts.F90
	gfortran -c -g ${FFLAGS} -o SURFEX/modd_ocean_csts.o                SURFEX/modd_ocean_csts.F90
	gfortran -c -g ${FFLAGS} -o SURFEX/modd_sltn.o                      SURFEX/modd_sltn.F90
	gfortran -c -g ${FFLAGS} -o SURFEX/modd_slt_surf.o                  SURFEX/modd_slt_surf.F90
	gfortran -c -g ${FFLAGS} -o SURFEX/modn_slt.o                       SURFEX/modn_slt.F90
	gfortran -c -g ${FFLAGS} -o SURFEX/modd_surfexn.o                   SURFEX/modd_surfexn.F90
	gfortran -c -g ${FFLAGS} -o SURFEX/init_slt.o                       SURFEX/init_slt.F90
	gfortran -c -g ${FFLAGS} -o SURFEX/coupling_sltn.o                  SURFEX/coupling_sltn.F90
	gfortran -c -g ${FFLAGS} -o SURFEX/mode_ssgf.o                      SURFEX/mode_ssgf.F90
	gfortran -c -g ${FFLAGS} -o driver_param_sea_spray_fluxes_offline.o driver_param_sea_spray_fluxes_offline.F90
	gfortran -v SURFEX/*.o *.o -o driver_param_sea_spray_fluxes_offline
	@echo
	@echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	@echo Successful compilation of
	@echo driver_param_sea_spray_fluxes_offline
	@echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	@echo
clean:
	rm -f ZSFSLT*
	rm -f *.mod *.o SURFEX/*.mod SURFEX/*.o
	rm -f driver_param_sea_spray_fluxes_offline
