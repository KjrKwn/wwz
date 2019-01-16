###### ######
use_gcc       = n
use_icc       = y
use_craygcc   = n
use_crayicc   = n
use_craycc    = n
debug         = n
optimize      = y
gprof         = n
use_omp       = n
use_mpi       = n
use_gsl       = n
use_hdf       = n
time_h        = n
###### ######
###### ######
use_xc = n
ifeq ($(ENV), xccfca)
use_xc = y
endif
ifeq ($(ENV), xc40)
use_xc = y
endif
# use_xc = n
###### ######

CC   = gcc
C++  = g++
FORT = gfortran
ifeq ($(use_icc), y)
CC   = icc
C++  = icpc
FORT = ifort
endif
ifeq ($(use_xc), y)
CC   = cc
C++  = CC
FORT = ftn
endif

ifeq ($(use_gcc), y)
CFLAGS += -lm
FFLAGS += -lm -fbounds-check -cpp -P
endif

ifeq ($(use_icc), y)
CFLAGS +=
FFLAGS += -check bounds -fpp
endif

ifeq ($(use_xc), y)
CFLAGS +=
FFLAGS += -fpp
endif

ifeq ($(debug), y)
CFLAGS += -DDEBUG
FFLAGS += -DDEBUG
ifeq ($(use_gcc), y)
CFLAGS += -g3
FFLAGS += -g3 -fbacktrace -dM
# FFLAGS += -Wall
endif
ifeq ($(use_craycc), y)
CFLAGS += -G0
FFLAGS += -G0
endif
ifeq ($(use_icc), y)
CFLAGS += -g -traceback
FFLAGS += -g -traceback -warn all
endif
endif


ifeq ($(optimize), y)
ifeq ($(use_gcc), y)
CFLAGS += -O3
FFLAGS += -O3
endif
ifeq ($(use_icc), y)
CFLAGS += -fast -parallel
FFLAGS += -fast -parallel
endif
ifeq ($(use_craygcc), y)
CFLAGS += -O3 -ftree-vectorize
FFLAGS += -O3 -ftree-vectorize
endif
ifeq ($(use_crayicc), y)
CFLAGS += -fast -parallel
FFLAGS += -fast -parallel
endif
ifeq ($(use_craycc), y)
CFLAGS += -O3 -h authothread
FFLAGS += -O3 -h authothread
endif
endif


ifeq ($(gprof), y)
CFLAGS += -pg
endif


ifeq ($(use_omp), y)
CFLAGS += -DOPENMP
FFLAGS += -DOPENMP
ifeq ($(use_gcc), y)
CFLAGS += -fopenmp
FFLAGS += -fopenmp
endif
ifeq ($(use_craycc), y)
CFLAGS += -h omp -h thread3
FFLAGS += -h omp -h thread3
endif
ifeq ($(use_crayicc), y)
CFLAGS += -qopenmp
FFLAGS += -qopenmp
endif
ifeq ($(use_icc), y)
CFLAGS += -qopenmp
FFLAGS += -qopenmp
endif
endif

ifeq ($(use_omp), n)
ifeq ($(use_craycc), y)
CFLAGS += -h noomp
FFLAGS += -h noomp
endif
endif


ifeq ($(use_mpi), y)
CFLAGS += -DMPI
FFLAGS += -DMPI
ifeq ($(ENV), local)
CC   = mpicc
C++  = mpicpp
FORT = mpifort
CLIBS  += -L/opt/local/lib/mpich-mp -lmpi
CFLAGS += -I/opt/local/include/mpich-mp
FLIBS  += -L/opt/local/lib/mpich-mp -lmpi
FFLAGS += -I/opt/local/include/mpich-mp
endif
endif


ifeq ($(use_gsl), y)
CFLAGS += -DGSL
ifeq ($(ENV), local)
GSL_LIBS += -L/opt/local/lib -lgsl -lgslcblas
GSL_CFLAGS += -I/opt/local/include
endif
endif


ifeq ($(use_hdf), y)
CFLAGS += -DHDF
ifeq ($(ENV), local)
CLIBS += -L/opt/local/lib -lhdf5 -lhdf5_hl
CFLAGS += -I/opt/local/include
endif
endif


ifeq ($(time_h), y)
CFLAGS +=-DWATCHTIME
endif


#########################################################################
############################         ####################################
############################ targets ####################################
############################         ####################################
#########################################################################
ETC = Makefile

libwwz.so: wwz_c.c $(ETC)
	$(CC) $(CLIBS) $(CFLAGS) -shared wwz_c.c -o libwwz.so

clean:
	# rm -f
	rm -f *~
	rm -f \#*\#

distclean:
	make clean
	rm -f *.o
	rm -f *.so
	rm -f *.pyc
	rm -f *.mod
	rm -f -r *.dSYM

allclean:
	make distclean
	rm -f run
