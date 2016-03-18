SUBDIRS := mbase \
                        mfbase \
                        mdata \
                        mfileio \
                        mhbase \
                        braw \
                        bbase \
                        bcalib \
                        bmcread \
                        bgeom \
                        bimpulse \
                        bjoint \
                        bsql \
                        bruninfo \
                        bfileio \
                        btrigger \
                        btelescope \
                        bhbase \
                        bhist \
                        bhvstime \
                        bhcalib  \
                        brawtreat \
                        bextractevent \
                        bdtcalib \
                        bstatistics \
                        bfilter \
                        bfqualify \
                        blasercalib

INCLUDES = -I. $(SUBDIRS:%=-I/home/local1/work/baikal/bars/BARS/tags/0.2.0/%) -I$(ROOTSYS)/include/root

ROOTLIBS   = $(shell root-config --libs)

LIBDIR = -L$(ROOTSYS)/lib/root -L. -L./lib -L${BW}

LIBS = 	-lmars \
	-lCore \
        -lCint \
	-lRIO \
        -lHist \
        -lGpad \
	-lTreePlayer \
        -lTree \
        -lRint \
        -lMatrix \
        -lPhysics \
        -lMathCore \
        -lThread \
	-lz \
       	-pthread \
       	-lm \
       	-ldl \
       	-rdynamic \
	-lstdc++ 

all: 
	g++ -v preprocessing.cc ${INCLUDES} ${LIBDIR} ${LIBS} 

