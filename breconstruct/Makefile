##################################################################
#
#   makefile
# 
#   for the MARS software
#
##################################################################
include ../Makefile.conf.general
include ../Makefile.conf.$(OSTYPE)

#------------------------------------------------------------------------------

#
# Handling name of the Root Dictionary Files
#
CINT  = BRec

#
#  connect the include files defined in the config.mk file
#
INCLUDES = -I. -I../mbase -I../mfileio -I../bgeom -I../bextractevent -I../bfilter -I../bfqualify -I../braw -I../bbase
# mgui (MCamEvent), mgeom(MGeomCam)

SRCFILES = 	BReconstruct.cc \
		BRecParameters.cc \
		BReconstructMuon.cc \
		BRecQualify.cc \
		BRecQuality.cc \
#		BReconstructLaser.cc \
############################################################

all: $(OBJS)

include ../Makefile.rules

mrproper:	clean rmbak
