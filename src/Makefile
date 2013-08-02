################################################################################
# Makefile fo DES ImagProc codes.
#
#  $Id: arcp 996 2008-05-12 20:31:51Z dadams $
#
#  $Rev:: 996                              $:  # Revision of last commit.
#  $LastChangedBy:: dadams                 $:  # Author of last commit. 
#  $LastChangedDate:: 2008-05-12 15:31:51 #$:  # Date of last commit.
#
#  The implit rule is for every *.o object file in the OBJ list, a cooresponding
#  *.c or *.cc file exists.  Executable targets are stil explixitly defined and
#  a complete list of working targets should be maintained in the EXE variable.
#
#  A future simplification would be to group list of main targets with the same
#  command and the same dependencies--reducing the number of explicit target definitions.
#
#  RAG added a line so that a new commit could supercede an erroneous checkin
################################################################################

################################################################################
# Install and include locations.
################################################################################
#
# working for backwards compat, locate libimsupport.a
#  Rule: if eups product use it.
#        if DES Home, use it (greg style stack)
#        n.b IMSUPPORT_DIR can be set with a setup -r .  from imsuport/trunk (or similar)
#
#

# $(info Using ${INCLUDE_FILE})
# $(warning DES_HOME is $(DES_HOME) AND ${DES_HOME})

ifeq ("$(IMSUPPORT_DIR)", "")
   LIBDETECT = $(wildcard $(DES_HOME)/lib/libimsupport.a))  
   ifneq ("(LIBDETECT)", "") 
       IMSUPPORT_DIR := $(DES_HOME)
   endif 
   $(info IM $(IMSUPPORT_DIR))
endif 
ifeq ("$(IMSUPPORT_DIR)", "")
   $(error imsupport not defined)
endif 
IMSUPPORT_LIB := ${IMSUPPORT_DIR}/lib/
IMSUPPORT_INCLUDE := ${IMSUPPORT_DIR}/include  
LIBIMSUPPORT:=${IMSUPPORT_LIB}/libimsupport.a

ifndef PREFIX
  ifdef DES_HOME
    PREFIX=${DES_HOME}
  else
    PREFIX=../..
 endif
endif

ifndef FITSIO_BASE
  ifdef DES_PREREQ
    FITSIO_BASE=${DES_PREREQ}
  else
    FITSIO_BASE=/usr
  endif
endif


ifndef NRECIPES_BASE
  ifdef DES_PREREQ
     NRECIPES_BASE=${DES_PREREQ}
  else
    NRECIPES_BASE=/usr/local
  endif
endif

ifndef FFTW_BASE
  ifdef DES_PREREQ
    FFTW_BASE=${DES_PREREQ}
  else
    FFTW_BASE = /usr/local
  endif
endif

FITSIO_LIB = $(FITSIO_BASE)/lib
FITSIO_INCLUDE = $(FITSIO_BASE)/cfitsio/include
CCFITS_INCLUDE = $(FITSIO_BASE)/include/CCfits
NRECIPES_LIB = $(NRECIPES_BASE)/lib
FFTW_LIB=$(FFTW_BASE)/lib
FFTW_INCLUDE=$(FFTW_BASE)/include



ifdef DES_PREREQ
  PSFEX_INCLUDE := ${DES_PREREQ}/include/psfex ${DES_PREREQ}/include/psfex/fits ${DES_PREREQ}/include/psfex/wcs ${DES_PREREQ}/include/fits ${DES_PREREQ}/include/wcs
  SWARP_INCLUDE := ${DES_PREREQ}/include/swarp
  PSFEX_LIB := ${DES_PREREQ}/lib/psfex
  SWARP_LIB := ${DES_PREREQ}/lib/swarp
  WCS_INCLUDE := ${DES_PREREQ}/include
  WCS_LIB := ${DES_PREREQ}/lib
  PLPLOT_LIB := ${DES_PREREQ}/lib
  PLPLOT_INCLUDE := ${DES_PREREQ}/include/plplot
else
  PSFEX_INCLUDE := ../../../terapix/trunk/psfex/src
  SWARP_INCLUDE := ../../../terapix/trunk/swarp/src
  PSFEX_LIB := ../../../terapix/trunk/psfex/src
  SWARP_LIB := ../../../terapix/trunk/swarp/src
endif

FITSIO :=       $(FITSIO_LIB)/libcfitsio.so
NR2 :=          $(NRECIPES_LIB)/libpress2ndkr.a
NR2D :=         $(NRECIPES_LIB)/libpress2ndkrD.a
FFTW := 	$(FFTW_LIB)/libfftw3f.a
FFTW_THREADS :=	$(FFTW_LIB)/libfftw3f_threads.a
#WCS_C :=        $(DES_PREREQ)/lib/psfex/wcs/libwcs_c.a
WCS_C :=        $(DES_PREREQ)/lib/wcs/libwcs_c.a
# TPIX_FITS :=    $(DES_PREREQ)/lib/psfex/fits/libfits.a
TPIX_FITS :=    $(DES_PREREQ)/lib/fits/libfits.a
#TPIX_FITS :=    $(PSFEX_LIB)/fits/libfits.a

BINDIR := 	$(PREFIX)/bin
LIBDIR := 	$(PREFIX)/lib

#INCLUDE :=    	$(FITSIO_INCLUDE) $(FFTW_INCLUDE) $(PSFEX_INCLUDE) $(SWARP_INCLUDE) ${WCS_INCLUDE}
#LIBS = $(FITSIO_LIB):$(NRECIPES_LIB):$(FFTW_LIB):$(PSFEX_LIB):$(SWARP_LIB):${WCS_LIB}

## Temp addition for WCSTOOLS library

DMWCS_INCLUDE := ${DES_PREREQ}/include
DMWCS_LIB := ${DES_PREREQ}/lib
DMWCS := ${DMWCS_LIB}/libwcs.a 

## Temp addition for PLPLOT library


INCLUDE :=    	$(CCFITS_INCLUDE) /usr/include/cfitsio $(FFTW_INCLUDE) $(PSFEX_INCLUDE) $(SWARP_INCLUDE) ${WCS_INCLUDE} ${DMWCS_INCLUDE} ${PLPLOT_INCLUDE} $(IMSUPPORT_INCLUDE) 
LIBS = $(FITSIO_LIB):$(NRECIPES_LIB):$(FFTW_LIB):$(PSFEX_LIB):$(SWARP_LIB):${WCS_LIB}:${DMWCS_LIB}:${PLPLOT_LIB}


# END Temp addition

# We could add all the headers to this list to force rebuilds 
# if a header changes
HEADERS := fitsio.h fft.h fitswcs.h imageproc.h

# Paths to search for dependecies:
vpath %.h $(INCLUDE)
vpath %.a $(LIBS)
vpath %.o $(LIBS)

################################################################################
# This is a default set of CC and CFLAGS settings.  they get overridden by anything
# passed in from a higher level make call
################################################################################
# gcc
CC :=		gcc -Wall
CPPC := 	g++
CPPFLAGS:=      -g -O3
CPPCOPTS :=	-c -Wall -funroll-loops -O3 -Wno-unused-variable 
FC :=		f77
CFLAGS := 	-O3 -g -Wall -Wno-unused-variable #-ftree-vectorize -ffast-math -msse2 -ftree-vectorizer-verbose=6  -DTRY2VECT
# Intel
#CC :=		icc
#CPPC := 	icpc
#CPPCOPTS :=	-O2 -w -g
#FC :=		ifort
#CFLAGS := 	-O2 -w -g

override CPPFLAGS += $(addprefix -I ,$(INCLUDE))
override CFLAGS += $(addprefix -I ,$(INCLUDE))


################################################################################
# Targets for all executables:
# Right now a target needs to be here for the executable to get built
# installed and cleaned up correctly.  The targets also need a build rule 
# defined.
# Separate lists can be defined if needed to define a alternate target, like 
# a set of related codes.
###############################################################################

#EXE_set1 = blah

EXE_REST =  ../bin/imcorrect			\
	../bin/DECam_crosstalk	                \
	../bin/mkillumcor	                \
	../bin/Mosaic2_crosstalk               \
	../bin/mkbiascor			\
	../bin/mkdarkcor			\
	../bin/mkflatcor			\
	../bin/mksupersky			\
        ../bin/mkbleedmask                     \
        ../bin/mkmask                          \
        ../bin/maskcosmics          

EXE = ${EXE_REST}
# Sample add another list:
# EXE = ${EXE_REST} ${EXE_set1}

# Sample add another target:
#.PHONY: somecodes
#somecodes: ${EXE_set1}

# Target to make all executables:
.PHONY: all
all: $(EXE) 


#############################################################
#############################################################
#          Begin object build rules section
#
# No object (*.o) targets should have to be explicitly defined
# unless there are special flags needed for that object.
#
#############################################################
#############################################################

# Special object targets all other obejcts will be built by
# the implicit pattern rule:
%.o : %.c $(HEADERS)
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
%.o : %.cc
	$(CPPC) $(CFLAGS) $(CPPFLAGS) $(CPPCOPTS) -c $^ -o $@




#############################################################
#############################################################
#          Begin Executable build rules section
#
#  Right now all rules are explicitly defined for all 
#  executable targets.  They have been grouped by dependencies
#  then external linkage.  At some point it may be worthwile 
#  to write a few rules that operate on a list of targets
#  that all share the same dependencies and linkage.
#
#############################################################
#############################################################

# "Solo" targets
# ../lib/imageutils.a 
../bin/DECam_crosstalk: DECam_crosstalk.o  ${LIBIMSUPPORT}
	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(FITSIO) -lm 
../bin/Mosaic2_crosstalk: Mosaic2_crosstalk.o  ${LIBIMSUPPORT}
	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(FITSIO) -lm 
../bin/imcorrect: imcorrect.o ${LIBIMSUPPORT}
	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(FITSIO) -lm
../bin/mkillumcor: mkillumcor.o ${LIBIMSUPPORT} 
	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(FITSIO) -lm
../bin/mksupersky:  mksupersky.o ${LIBIMSUPPORT}   
	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(FITSIO) -lm
../bin/mkbiascor:  mkbiascor.o ${LIBIMSUPPORT}   
	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(FITSIO) -lm
../bin/mkdarkcor:  mkdarkcor.o ${LIBIMSUPPORT}   
	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(FITSIO) -lm
../bin/mkflatcor:  mkflatcor.o ${LIBIMSUPPORT}   
	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(FITSIO) -lm
../bin/mkbpm:  mkbpm.o ${LIBIMSUPPORT}  
	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(FITSIO) -lm
../bin/mkbleedmask: MakeBleedMask.o BleedTrail.o ${LIBIMSUPPORT} fitswcs.o   
	$(CPPC) -o $@ $^ -L$(DES_PREREQ)/lib $(WCS_C) $(TPIX_FITS) $(FITSIO) -lCCfits -lm
../bin/mkmask: mkmask.o mkstarmask.o fitswcs.o ${LIBIMSUPPORT} 
	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(WCS_C) $(TPIX_FITS) $(FITSIO) -lm
../bin/maskcosmics: maskcosmics.o mkstarmask.o fitswcs.o ${LIBIMSUPPORT} 
	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(WCS_C) $(TPIX_FITS) $(FITSIO) -lm


################################################################################
# Clean-up rules.
################################################################################
tidy:
	rm -rf *.o *~ \#*\#

clean: tidy
	for f in $(EXE_REST) ; do rm -f $$f ; done


install:
	false #nothing to install here.
