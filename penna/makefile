.PHONY: all clean test

HPPFILES = 	src/genome.hpp 		\
			src/animal.hpp 		\
			src/population.hpp  \
			src/fishing_pop.hpp

CXXFLAGS  += -Isrc
LDFLAGS	  += -Lsrc
LDLIBS 	  += -lpenna
MAKEFLAGS += --no-print-directory

all: pop.pdf

libpenna.a: ${HPPFILES}
	${MAKE} -Csrc all

test: libpenna.a
	${MAKE} -Ctest all
 
main.o: main.cpp src/fishing_pop.hpp test
	${CXX} ${CXXFLAGS} -c $< -o $@ 

main.exe: main.o libpenna.a
	${CXX} ${CXXFLAGS} $< -o $@ ${LDFLAGS} ${LDLIBS}

data.txt: main.exe
	./$<

pop.pdf: plot.py data.txt
	python $<

clean:
	${MAKE} -Csrc clean
	${RM} *.{o,a,exe,pdf} data.txt