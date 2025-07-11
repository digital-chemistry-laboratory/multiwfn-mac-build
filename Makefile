SIMD = -msse3
DISDIAG = -diag-disable 8290,8291,6371,10316,6178,6916,7416,5268,7342,7373,5120,5144,5082,5112,2554,5183,6182,7352
OPT = -O2 -qopenmp -qopenmp-link=static -threads -qopt-matmul $(SIMD) $(DISDIAG) -fpscomp logicals -fpp -mkl -static-intel -DINTEL_MKL -stand f08
OPT1 = -O1 -qopenmp -qopenmp-link=static -threads $(SIMD) $(DISDIAG) -fpscomp logicals -fpp -mkl -static-intel -DINTEL_MKL -stand f08
#Options in the next line is for debugging purpose
#OPT = -O0 -qopenmp -qopenmp-link=static -threads $(DISDIAG) -fpscomp logicals -fpp -mkl -static-intel -DINTEL_MKL -stand f08 -debug all -g -traceback -check all -fstack-protector

LIB_base = 
LIB_GUI = $(LIB_base) ./dislin_d-11.0.a -lXm -lXt -lX11 -lGL
LIB_noGUI = $(LIB_base)
INCLUDE = -I./ -I./ext
FC = ifort
CC = gcc
EXE = Multiwfn
EXE_noGUI = Multiwfn_noGUI
LIBRETAPATH = ./libreta_hybrid

objects = define.o util.o plot.o Bspline.o sym.o libreta.o function.o GUI.o sub.o integral.o Lebedev-Laikov.o \
DFTxclib.o edflib.o fparser.o fileIO.o spectrum.o DOS.o Multiwfn.o 0123dim.o LSB.o \
population.o frj.o orbcomp.o bondorder.o topology.o excittrans.o otherfunc.o \
otherfunc2.o otherfunc3.o O1.o surfana.o procgriddata.o AdNDP.o fuzzy.o CDA.o basin.o \
orbloc.o visweak.o EDA.o CDFT.o ETS_NOCV.o atmraddens.o NAONBO.o grid.o PBC.o hyper_polar.o deloc_aromat.o cp2kmate.o\
minpack.o blockhrr_012345.o ean.o hrr_012345.o eanvrr_012345.o boysfunc.o naiveeri.o ryspoly.o 2F2.f90.o

objects_noGUI = noGUI/dislin_d_empty.o noGUI/mouse_rotate_empty.o #Dummy subroutines for noGUI version

ifeq ($(WITH_FD),1)
  objects += 2F2.c.o
  ifeq ($(OS),Ubuntu) # for Ubuntu
    LIB_base += -lflint -lflint-arb
  else ifeq ($(OS),RHEL) # for Fedora, CentOS, RHEL
    INCLUDE += -I/usr/include/arb
    LIB_base += -lflint -larb
  endif
else
  objects += no2F2.c.o
endif

default: $(objects)
	$(MAKE) noGUI
	$(MAKE) GUI
	@echo " ------------------------------------------------------ "
	@echo "          Multiwfn has been successfully built!"
	@echo " ------------------------------------------------------ "

GUI: $(objects) mouse_rotate.o xlib.o
	$(FC) $(OPT) $(objects) mouse_rotate.o xlib.o $(LIB_GUI) -o $(EXE)

noGUI: $(objects) $(objects_noGUI)
	$(FC) $(OPT) $(objects) $(objects_noGUI) $(LIB_noGUI) -o $(EXE_noGUI)

clean:
	rm -f $(EXE) *.o *.mod

#Only clean Multiwfn files, compiled libreta files are not affected
cleanmultiwfn:
	mkdir tmplib
	mv libreta.o ean.o hrr_012345.o blockhrr_012345.o eanvrr_012345.o boysfunc.o \
	libreta.mod hrr.mod blockhrr.mod ean.mod eanvrr.mod boysfunc.mod tmplib
	rm -f $(EXE) *.o *.mod
	mv tmplib/* ./
	rm -r tmplib

#Only clean libreta files, Multiwfn libreta files are not affected
cleanlibreta:
	rm -f $(EXE) libreta.o ean.o hrr_012345.o blockhrr_012345.o eanvrr_012345.o \
	boysfunc.o libreta.mod hrr.mod blockhrr.mod ean.mod eanvrr.mod boysfunc.mod


#Define modules that used by other files

dislin.mod : dislin_d.f90
	$(FC) $(OPT) -c dislin_d.f90

define.o : define.f90 dislin.mod
	$(FC) $(OPT) -c define.f90

Bspline.o : Bspline.f90
	$(FC) $(OPT) -c Bspline.f90

util.o : util.f90 define.o
	$(FC) $(OPT) -c util.f90

function.o : function.f90 define.o util.o Bspline.o libreta.o 2F2.f90.o
	$(FC) $(OPT) -c function.f90

plot.o : plot.f90 function.o define.o util.o
	$(FC) $(OPT) -c plot.f90

GUI.o : GUI.f90 define.o plot.o function.o mouse_rotate.o
	$(FC) $(OPT) -c GUI.f90

mouse_rotate.o : mouse_rotate.f90 xlib.o define.o plot.o
	$(FC) $(OPT) -c mouse_rotate.f90

noGUI/mouse_rotate_empty.o : noGUI/mouse_rotate_empty.f90
	$(FC) $(OPT) -c noGUI/mouse_rotate_empty.f90 -o noGUI/mouse_rotate_empty.o

2F2.f90.o : ext/2F2.f90 util.o Bspline.o
	$(FC) $(OPT) -c ext/2F2.f90 -o 2F2.f90.o

modules = define.o util.o function.o plot.o GUI.o libreta.o 2F2.f90.o


#Library or adpated third-part codes

DFTxclib.o : DFTxclib.F
	$(FC) $(OPT) -c DFTxclib.F

Lebedev-Laikov.o : Lebedev-Laikov.F
	$(FC) $(OPT) -c Lebedev-Laikov.F

sym.o : sym.F
	$(FC) $(OPT) -c sym.F

edflib.o : edflib.f90
	$(FC) $(OPT) -c edflib.f90

atmraddens.o : atmraddens.f90
	$(FC) $(OPT) -c atmraddens.f90

minpack.o : minpack.f90
	$(FC) $(OPT) -c minpack.f90
	
fparser.o : fparser.f90
	$(FC) $(OPT) -c fparser.f90
	
2F2.c.o : ext/2F2.c
	$(CC) $(INCLUDE) -c ext/2F2.c -o 2F2.c.o

no2F2.c.o : ext/no2F2.c
	$(CC) -c ext/no2F2.c -o no2F2.c.o

noGUI/dislin_d_empty.o : noGUI/dislin_d_empty.f90
	$(FC) $(OPT) -c noGUI/dislin_d_empty.f90 -o noGUI/dislin_d_empty.o -diag-disable 6178,6843

#Others

sub.o : sub.f90 $(modules)
	$(FC) $(OPT) -c sub.f90

integral.o : integral.f90 $(modules)
	$(FC) $(OPT) -c integral.f90

fileIO.o : fileIO.f90 $(modules)
	$(FC) $(OPT) -c fileIO.f90

spectrum.o : spectrum.f90 $(modules)
	$(FC) $(OPT) -c spectrum.f90

DOS.o : DOS.f90 $(modules)
	$(FC) $(OPT) -c DOS.f90

Multiwfn.o : Multiwfn.f90 $(modules)
	$(FC) $(OPT) -c Multiwfn.f90

0123dim.o : 0123dim.f90 $(modules)
	$(FC) $(OPT) -c 0123dim.f90

LSB.o : LSB.f90 $(modules)
	$(FC) $(OPT) -c LSB.f90

population.o : population.f90 $(modules)
	$(FC) $(OPT) -c population.f90

frj.o : ext/frj.f90 $(modules)
	$(FC) $(OPT) -c ext/frj.f90 -o frj.o

orbcomp.o : orbcomp.f90 $(modules)
	$(FC) $(OPT) -c orbcomp.f90

bondorder.o : bondorder.f90 $(modules)
	$(FC) $(OPT) -c bondorder.f90

topology.o : topology.f90 $(modules)
	$(FC) $(OPT) -c topology.f90

excittrans.o : excittrans.f90 $(modules)
	$(FC) $(OPT) -c excittrans.f90

otherfunc.o : otherfunc.f90 $(modules)
	$(FC) $(OPT) -c otherfunc.f90

otherfunc2.o : otherfunc2.f90 $(modules)
	$(FC) $(OPT) -c otherfunc2.f90

otherfunc3.o : otherfunc3.f90 $(modules)
	$(FC) $(OPT) -c otherfunc3.f90
	
O1.o : O1.f90 $(modules)
	$(FC) $(OPT1) -c O1.f90

surfana.o : surfana.f90 $(modules)
	$(FC) $(OPT) -c surfana.f90

procgriddata.o : procgriddata.f90 $(modules)
	$(FC) $(OPT) -c procgriddata.f90

AdNDP.o : AdNDP.f90 $(modules)
	$(FC) $(OPT) -c AdNDP.f90

fuzzy.o : fuzzy.f90 $(modules)
	$(FC) $(OPT) -c fuzzy.f90

CDA.o : CDA.f90 $(modules)
	$(FC) $(OPT) -c CDA.f90

basin.o : basin.f90 $(modules)
	$(FC) $(OPT) -c basin.f90

orbloc.o : orbloc.f90 $(modules)
	$(FC) $(OPT) -c orbloc.f90

visweak.o : visweak.f90 $(modules)
	$(FC) $(OPT) -c visweak.f90

EDA.o : EDA.f90 $(modules)
	$(FC) $(OPT) -c EDA.f90

CDFT.o : CDFT.f90 $(modules)
	$(FC) $(OPT) -c CDFT.f90

ETS_NOCV.o : ETS_NOCV.f90 $(modules)
	$(FC) $(OPT) -c ETS_NOCV.f90

NAONBO.o : NAONBO.f90 $(modules)
	$(FC) $(OPT) -c NAONBO.f90

grid.o : grid.f90 $(modules)
	$(FC) $(OPT) -c grid.f90

PBC.o : PBC.f90 $(modules)
	$(FC) $(OPT) -c PBC.f90

hyper_polar.o : hyper_polar.f90 $(modules)
	$(FC) $(OPT) -c hyper_polar.f90
	
deloc_aromat.o : deloc_aromat.f90 $(modules)
	$(FC) $(OPT) -c deloc_aromat.f90
	
cp2kmate.o : cp2kmate.f90 $(modules)
	$(FC) $(OPT) -c cp2kmate.f90


# Interfaces of libreta-ESP to Multiwfn

libreta.o: ${LIBRETAPATH}/libreta.f90 hrr_012345.o blockhrr_012345.o ean.o eanvrr_012345.o boysfunc.o
	$(FC) $(OPT) -c ${LIBRETAPATH}/libreta.f90


# Pure libreta files for ESP

hrr_012345.o: ${LIBRETAPATH}/hrr_012345.f90
	$(FC) $(OPT) -diag-disable 6843 $(SIMD) -c ${LIBRETAPATH}/hrr_012345.f90

blockhrr_012345.o: ${LIBRETAPATH}/blockhrr_012345.f90
	$(FC) -O1 -diag-disable 6843 $(SIMD) -c ${LIBRETAPATH}/blockhrr_012345.f90

ean.o: ${LIBRETAPATH}/ean.f90 hrr_012345.o eanvrr_012345.o boysfunc.o ${LIBRETAPATH}/ean_data1.h ${LIBRETAPATH}/ean_data2.h
	$(FC) $(OPT) -c ${LIBRETAPATH}/ean.f90

eanvrr_012345.o: ${LIBRETAPATH}/eanvrr_012345.f90 boysfunc.o
	$(FC) $(OPT) -c ${LIBRETAPATH}/eanvrr_012345.f90

boysfunc.o: ${LIBRETAPATH}/boysfunc.f90 ${LIBRETAPATH}/boysfunc_data1.h ${LIBRETAPATH}/boysfunc_data2.h
	$(FC) $(OPT) -c ${LIBRETAPATH}/boysfunc.f90


# libreta-ERI

naiveeri.o: ${LIBRETAPATH}/naiveeri.f90 ryspoly.o
	$(FC) $(OPT) -c ${LIBRETAPATH}/naiveeri.f90
	
ryspoly.o: ${LIBRETAPATH}/ryspoly.f90
	$(FC) $(OPT) -c ${LIBRETAPATH}/ryspoly.f90


# Fortran-xlib interface
xlib.o: ext/xlib.f90
	$(FC) $(OPT) -fpscomp logicals -c ext/xlib.f90

