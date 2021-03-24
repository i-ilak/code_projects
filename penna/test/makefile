.PHONY: all clean testdone

SOURCEDIR  = ./../src/
SOURCELIB  = ${SOURCEDIR}libpenna.a
SOURCEFILE = ${SOURCEDIR}genome.hpp
LDLIBS 	  += -lpenna
LDFLAGS   += -L${SOURCEDIR}
HEADERSRC  = -I${SOURCEDIR}

all: testdone

copy_ctor.o: copy_ctor.cpp ${SOURCEFILE}
	${CXX} ${CXXFLAGS} -c $< -o $@ ${HEADERSRC}

genome_default_ctor.o: genome_default_ctor.cpp ${SOURCEFILE}
	${CXX} ${CXXFLAGS} -c $< -o $@ ${HEADERSRC}

test1.exe: copy_ctor.o  ${SOURCELIB}
	${CXX} ${CXXFLAGS} -o $@ $< ${LDLIBS} ${LDFLAGS}

test2.exe: genome_default_ctor.o ${SOURCELIB}
	${CXX} ${CXXFLAGS} -o $@ $< ${LDLIBS} ${LDFLAGS}

testdone: test1.exe test2.exe
	./$<
	./test2.exe

clean:
	${RM} *.{o,exe,c}