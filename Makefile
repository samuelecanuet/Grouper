ROOT_CFLAGS     = $(shell root-config --cflags)
ROOT_LIBS       = $(shell root-config --libs)

GSL_CFLAGS      = $(shell gsl-config --cflags)
GSL_LIBS        = $(shell gsl-config --libs)

GTOOLS_L = $(shell gtools_lib)

C++_CFLAGS      = $(shell -std=c++17)


SignalDict.cxx: Signal.h LinkDef.h
	rootcling -f $@ -c $^ 

CC        = g++
CFLAGS    = ${FASTERAC_CFLAGS} ${ROOT_CFLAGS} ${GSL_CFLAGS}
LIBS      = ${FASTERAC_LIBS}   ${ROOT_LIBS} ${GSL_LIBS} ${GTOOLS_L}

SRCEXE    = $(shell ls *.c++)
EXE       = $(SRCEXE:.c++=)
SOURCE    = $(shell ls *.hh)

all : $(EXE) $(SOURCE) 

$(EXE): $(SOURCE) $(SRCEXE)
	${CC} $@.c++ $@.hh Detectors.hh SignalDict.cxx -o $@ ${CFLAGS} ${LIBS} -lboost_filesystem -lboost_system

clean :
	rm -f *.o
	rm -f $(EXE)
	rm -f SignalDict.cxx
	rm -f SignalDict_rdict.pcm


