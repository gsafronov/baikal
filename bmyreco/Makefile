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
CINT  = BMyRec

#
#  connect the include files defined in the config.mk file
#
INCLUDES = -I. -I../mbase -I../mfileio -I../bgeom -I../bextractevent -I../bfilter -I../bfqualify -I../braw -I../bbase -I../bimpulse -I../bmcread -I../breconstruct
# mgui (MCamEvent), mgeom(MGeomCam)

SRCFILES = 	BMyReco.cc \
		BMyRecoMC.cc \
		BMyRecoBEvt.cc \
		BMyRecoReco.cc
############################################################

all: $(OBJS)

include ../Makefile.rules

mrproper:	clean rmbak
