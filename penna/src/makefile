.PHONY: all clean ofiles

AR 		 = ar
ARFLAGS  = rcs
HPPFILES = genome.hpp animal.hpp population.hpp fishing_pop.hpp
SRCFILES = ${HPPFILES:.hpp=.cpp}
OFILES   = ${HPPFILES:.hpp=.o}

all: libpenna.a

ofiles: ${HPPFILES}
	${CXX} ${CXXFLAGS} -c ${SRCFILES}

libpenna.a: ${OFILES}
	${AR} ${ARFLAGS} $@ ${OFILES}

clean: 
	${RM} *.{o,a,hpp.gch}