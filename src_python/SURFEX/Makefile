##########################################################
#
# To compile the code :
# make
#
# To clean the compile directory :
# make clean
#
##########################################################

FFLAGS=-O0 -fimplicit-none -Wall
CPPFLAGS=-I${PWD}

default:
	gfortran -c -g ${FFLAGS} -o modd_csts.o       modd_csts.F90
	gfortran -c -g ${FFLAGS} -o modd_ocean_csts.o modd_ocean_csts.F90
	f2py -m mode_ssgf -c mode_ssgf.F90 modd_csts.F90 modd_ocean_csts.F90 
	@echo
	@echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	@echo Successful compilation of
	@echo mode_ssgf
	@echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	@echo
clean:
	rm -f *.mod *.o *.so
