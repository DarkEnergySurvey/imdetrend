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

# Simplified to make it EUPS compliant -- FM

# IMSUPPORT_DIR is defined via EUPS
ifndef IMSUPPORT_DIR
   $(error imsupport not defined)
endif 

IMSUPPORT_LIB     := ${IMSUPPORT_DIR}/lib/
IMSUPPORT_INCLUDE := ${IMSUPPORT_DIR}/include  
LIBIMSUPPORT      := ${IMSUPPORT_LIB}/libimsupport.a

ifndef PREFIX
  PREFIX=../..
endif

# -------------------------------------------------
# Path setup for eups 1.2.30 package.
# Define Flavor Shared lib extension -- FM
FLAVOR = $(shell uname -s) # this give either Darwin/Linux
$(info Discovered FLAVOR: $(FLAVOR))
ifneq (,$(findstring Darwin, ${FLAVOR}))
   SHLIB_SUFFIX = .dylib
endif 
ifneq (,$(findstring Linux, ${FLAVOR}))
   SHLIB_SUFFIX = .so
endif 
INCLUDE   := $(CFITSIO_SHARED_DIR)/include $(CCFITS_DIR)/include $(CCFITS_DIR)/include/CCfits $(FFTW_DIR)/include $(IMSUPPORT_DIR)/include $(SWARP_DIR)/include
LIBS       = $(NRECIPES_DIR)/lib $(FFTW_DIR)/lib $(SWARP_DIR)/lib $(CFITSIO_SHARED_DIR)/lib
TPIX_FITS := $(SWARP_DIR)/lib/fits/libfits.a
FITSIO    := -lcfitsio
CCFITS    := $(CCFITS_DIR)/lib/libCCfits$(SHLIB_SUFFIX) 
$(info Setting up CCFITS: $(CCFITS) for $(FLAVOR) )
WCS_C     := $(SWARP_DIR)/lib/wcs/libwcs_c.a
NR2       := $(NRECIPES_DIR)/lib/libpress2ndkr.a
# END of EUPS addition
# --------------------------------------------------- FM

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
CPPC := 	$(CXX) -std=gnu++98
CPPFLAGS:=      -g -O3
CPPCOPTS :=	-c -Wall -funroll-loops -O3 -Wno-unused-variable 
FC :=		f77
CFLAGS := 	-O3 -g -Wall -Wno-unused-variable #-ftree-vectorize -ffast-math -msse2 -ftree-vectorizer-verbose=6  -DTRY2VECT

override CPPFLAGS += $(addprefix -I ,$(INCLUDE))
override CFLAGS += $(addprefix -I ,$(INCLUDE))


override CPPFLAGS += $(addprefix -L ,$(LIBS))
override CFLAGS += $(addprefix -L ,$(LIBS))
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
        ../bin/mkbpm

#
# ADW 2016-06-09: removed maskcosmics (and the need for BADPIX_LOW)
#
#        ../bin/maskcosmics          		\

#
# RAG removed mkmask (and therefore the need to compile mkstarmask
#
#        ../bin/mkmask                          \

EXE = ${EXE_REST}

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
	$(CPPC) $(CFLAGS) -o $@ $^ -L$(DES_PREREQ)/lib $(WCS_C) $(TPIX_FITS) $(FITSIO) $(CCFITS) -lm

# ADW 2016-06-09: removed
#../bin/maskcosmics: maskcosmics.o fitswcs.o ${LIBIMSUPPORT} 
#	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(WCS_C) $(TPIX_FITS) $(FITSIO) -lm

# RAG removed 
#../bin/mkmask: mkmask.o mkstarmask.o fitswcs.o ${LIBIMSUPPORT} 
#	$(CC) $(CFLAGS) -o $@ $^  $(NR2) $(WCS_C) $(TPIX_FITS) $(FITSIO) -lm


################################################################################
# Clean-up rules.
################################################################################
tidy:
	rm -rf *.o *~ \#*\#

clean: tidy
	for f in $(EXE_REST) ; do rm -f $$f ; done

install:
	false #nothing to install here.
