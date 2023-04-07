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
CPPFLAGS=-I${PWD}

default:
	gfortran -c -g ${FFLAGS} -o modd_csts.o modd_csts.F90
	gfortran -c -g ${FFLAGS} -o modd_sltn.o modd_sltn.F90
	gfortran -c -g ${FFLAGS} -o modd_slt_surf.o modd_slt_surf.F90
	gfortran -c -g ${FFLAGS} -o modd_surfexn.o modd_surfexn.F90
	gfortran -c -g ${FFLAGS} -o init_slt.o init_slt.F90
	gfortran -c -g ${FFLAGS} -o coupling_sltn.o coupling_sltn.F90
	gfortran -c -g ${FFLAGS} -o driver_param_sea_spray_fluxes_offline.o driver_param_sea_spray_fluxes_offline.F90
	gfortran -v *.o -o driver_param_sea_spray_fluxes_offline	
clean:
	rm -f *.mod *.o *_.f *_.F90
	rm -f driver_param_sea_spray_fluxes_offline
